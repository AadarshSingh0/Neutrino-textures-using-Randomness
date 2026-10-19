(* TheorySpaceNeutrinos 1.0. Real couplings only. All dimensional inputs in TeV.
   No experimental cuts, octant folding, or implicit changes of other parameters. *)
BeginPackage["TheorySpaceNeutrinos`"];
TSParameters::usage = "TSParameters[regime, geometry, nature, channel] returns an editable parameter Association. Channel is Hamiltonian or Boundary; Majorana is an alias for Boundary for Majorana neutrinos.";
TSGenerateBlocks::usage = "TSGenerateBlocks[p] returns six reproducible, exact symmetric tower blocks.";
TSMass::usage = "TSMass[p] generates six blocks and returns masses, active mixing and diagnostics.";
TSMassFromBlocks::usage = "TSMassFromBlocks[blocks,p] solves an explicitly supplied six-block realization, bypassing random generation.";
TSScan::usage = "TSScan[p,n,\"Parallel\"->False] keeps all n samples, with seeds p[\"Seed\"] through p[\"Seed\"]+n-1.";
TSPrecisionCheck::usage = "TSPrecisionCheck[p,extraDigits:20] compares two precisions on identical exact blocks.";
TSMixingAngles::usage = "TSMixingAngles[u] returns {theta12,theta13,theta23} in degrees from a 3 by 3 active block. These are raw angle diagnostics if u is nonunitary.";
TSPetersenAdjacency::usage = "TSPetersenAdjacency[k] returns the paper's generalized Petersen labeling, with inner vertices first. k must be divisible by 4 and at least 8.";
TSTable::usage = "TSTable[results] returns a compact Dataset; full results retain mixing matrices and diagnostics.";
Begin["`Private`"];

$packageFile = ExpandFileName[$InputFileName];
$version = "1.0";
$blockMap = {{1,2,3},{2,4,5},{3,5,6}};
$weakH = {{7/50,1,3/10},{1,1/2,7/20},{2/5,7/10,1}};
$weakY = {{1,1/2,2/5},{1/2,1,3/10},{1/2,9/10,1}};
$strongFlavor = {{1,2/5,3/5},{3/10,1,4/5},{9/10,3/10,1}};
fail[tag_, text_] := Throw[Failure[tag,<|"MessageTemplate"->text|>], tsFailure];
need[test_,tag_,text_] := If[!TrueQ[test],fail[tag,text]];
maxAbs[x_] := Max[Abs[Flatten[x]]];
exact[x_] := Rationalize[x,0];
realQ[x_] := NumericQ[x] && TrueQ[Im[x]==0];
matrix3Q[x_] := MatrixQ[x,realQ] && Dimensions[x]=={3,3};

TSParameters[regime_String:"WeakSite", geometry_String:"Local",
    nature_String:"Majorana", channel_String:"Hamiltonian"] := Catch[Module[
  {p, f, y, boundary, k},
  need[MemberQ[{"WeakSite","StrongSite","WeakHopping","StrongHopping","Custom"},regime],"Regime","Unknown regime."];
  need[MemberQ[{"Local","NonLocal","Petersen"},geometry],"Geometry","Unknown geometry."];
  need[MemberQ[{"Dirac","Majorana"},nature],"Nature","Nature must be Dirac or Majorana."];
  need[MemberQ[{"Hamiltonian","Boundary"},channel] || (nature=="Majorana" && channel=="Majorana"),"Channel","Use Hamiltonian or Boundary (Majorana is an alias for Boundary in the Majorana case)."];
  boundary=channel!="Hamiltonian";
  f=If[regime=="StrongSite",$strongFlavor,$weakH];
  y=If[regime=="StrongSite",$strongFlavor,$weakY];
  k=If[regime=="WeakSite",If[geometry=="Local",9,16],8];
  p=<|"Nature"->nature,"Geometry"->geometry,"Preset"->regime,"Channel"->If[boundary,"Boundary","Hamiltonian"],
    "Sites"->k,"SiteMean"->10,"SiteHalfWidth"->1/5,
    "HoppingMean"->1/5,"HoppingHalfWidth"->0,"DecayBase"->3,
    "DistanceNormalization"->"NearestNeighbor","EdgeLaw"->"Uniform",
    "HamiltonianFlavor"->If[boundary,IdentityMatrix[3],f],
    "ActiveCoupling"->173/1000,"ActiveFlavor"->IdentityMatrix[3],
    "DiracRightCoupling"->173/1000,"DiracRightFlavor"->If[boundary,y,IdentityMatrix[3]],
    "PsiCoupling"->1/5,"PsiFlavor"->IdentityMatrix[3],
    "MajoranaScale"->10,"MajoranaFlavor"->If[boundary,y,IdentityMatrix[3]],
    "SymmetrizeMajorana"->True,"WorkingPrecision"->80,"RandomDigits"->24,
    "Seed"->20261003,"Method"->"Rotated","MassOrdering"->"Ascending",
    "KeepMatrices"->False,"DegeneracyTolerance"->10^-12,"ResidualTolerance"->10^-8|>;
  Switch[regime,
    "StrongSite",p=Join[p,<|"SiteMean"->0,"SiteHalfWidth"->10,"MajoranaScale"->5|>],
    "WeakHopping",p=Join[p,<|"SiteHalfWidth"->0,"HoppingMean"->0,"HoppingHalfWidth"->1/10,
      "DecayBase"->2,"DistanceNormalization"->"Numerator","EdgeLaw"->"LegacyAverage","PsiCoupling"->1/10|>],
    "StrongHopping",p=Join[p,<|"SiteMean"->Missing["SetExplicitly"],"SiteHalfWidth"->0,
      "HoppingMean"->0,"HoppingHalfWidth"->Missing["SetExplicitly"],"DecayBase"->2,
      "DistanceNormalization"->"Numerator","EdgeLaw"->"Uniform",
      "MajoranaScale"->Missing["SetExplicitly"],"PsiCoupling"->Missing["SetExplicitly"]|>],
    "Custom",Null];
  p],tsFailure];

validate[p_Association] := Module[{keys,missing,unknown,scalarKeys,matrixKeys},
  keys=Keys[TSParameters[]]; missing=Complement[keys,Keys[p]];unknown=Complement[Keys[p],keys];
  need[missing=={},"Parameters","Missing keys: "<>ToString[missing,InputForm]];
  need[unknown=={},"Parameters","Unknown keys (possible typo): "<>ToString[unknown,InputForm]];
  need[MemberQ[{"Dirac","Majorana"},p["Nature"]],"Nature","Use Dirac or Majorana."];
  need[MemberQ[{"Local","NonLocal","Petersen"},p["Geometry"]],"Geometry","Unknown geometry."];
  need[IntegerQ[p["Sites"]] && p["Sites"]>=2,"Sites","Sites must be an integer at least 2."];
  If[p["Geometry"]=="Petersen",need[p["Sites"]>=8 && Mod[p["Sites"],4]==0,"Sites","Petersen requires Sites divisible by 4 and at least 8."]];
  scalarKeys={"SiteMean","SiteHalfWidth","HoppingMean","HoppingHalfWidth","DecayBase",
    "ActiveCoupling","DiracRightCoupling","PsiCoupling","MajoranaScale","DegeneracyTolerance","ResidualTolerance"};
  Do[need[realQ[p[key]],"RealInput","Set a real numeric value for "<>key<>". Complex couplings are not supported in this release."],{key,scalarKeys}];
  need[p["SiteHalfWidth"]>=0 && p["HoppingHalfWidth"]>=0,"Widths","Half widths must be nonnegative."];
  need[p["DecayBase"]>0,"DecayBase","DecayBase must be positive."];
  need[p["DegeneracyTolerance"]>=0 && p["ResidualTolerance"]>0,"Tolerance","Invalid diagnostic tolerance."];
  matrixKeys={"HamiltonianFlavor","ActiveFlavor","DiracRightFlavor","PsiFlavor","MajoranaFlavor"};
  Do[need[matrix3Q[p[key]],"FlavorMatrix",key<>" must be a real 3 by 3 matrix."],{key,matrixKeys}];
  need[IntegerQ[p["WorkingPrecision"]] && p["WorkingPrecision"]>=30,"Precision","WorkingPrecision must be an integer at least 30."];
  need[IntegerQ[p["RandomDigits"]] && p["RandomDigits"]>=1,"RandomDigits","RandomDigits must be a positive integer."];
  need[IntegerQ[p["Seed"]],"Seed","Seed must be an integer."];
  Do[need[MemberQ[{True,False},p[key]],"Boolean",key<>" must be True or False."],{key,{"KeepMatrices","SymmetrizeMajorana"}}];
  need[MemberQ[{"NearestNeighbor","Numerator"},p["DistanceNormalization"]],"Normalization","Use NearestNeighbor or Numerator."];
  need[MemberQ[{"Uniform","LegacyAverage"},p["EdgeLaw"]],"EdgeLaw","Use Uniform or LegacyAverage."];
  need[MemberQ[{"Direct","Rotated"},p["Method"]],"Method","Use Direct or Rotated."];
  need[MemberQ[{"Ascending","Inverted"},p["MassOrdering"]],"MassOrdering","Use Ascending or Inverted."];
  If[p["Nature"]=="Majorana" && !p["SymmetrizeMajorana"],
    need[TrueQ[exact[p["MajoranaFlavor"]]==Transpose[exact[p["MajoranaFlavor"]]]],"MajoranaSymmetry","MajoranaFlavor must be symmetric, or explicitly enable SymmetrizeMajorana."]];
  p
];

TSPetersenAdjacency[k_Integer] := Catch[Module[{a,n,pairs},
  need[k>=8 && Mod[k,4]==0,"Sites","Petersen requires k divisible by 4 and at least 8."];
  n=k/2;a=ConstantArray[0,{k,k}];
  Do[pairs={{n+i,n+Mod[i,n]+1},{i,n+i},{i,Mod[i-1+n/2,n]+1}};
    Do[a[[pair[[1]],pair[[2]]]]=1;a[[pair[[2]],pair[[1]]]]=1,{pair,pairs}],{i,n}];a],tsFailure];

generate[p_] := Module[{k=p["Sites"],den=10^p["RandomDigits"],sites,r1,r2,adj,blocks,q,d},
  (* Draw all potential edges before geometry masking, to permit matched seeds. *)
  BlockRandom[SeedRandom[p["Seed"],Method->"MersenneTwister"];
    sites=RandomInteger[{-den,den},{6,k}]/den;
    r1=RandomInteger[{-den,den},{6,k,k}]/den;
    r2=RandomInteger[{-den,den},{6,k,k}]/den];
  adj=If[p["Geometry"]=="Petersen",TSPetersenAdjacency[k],ConstantArray[1,{k,k}]];
  blocks=Table[DiagonalMatrix[exact[p["SiteMean"]]+exact[p["SiteHalfWidth"]]*sites[[b]]],{b,6}];
  Do[d=j-i;
    q=exact[p["HoppingMean"]]+exact[p["HoppingHalfWidth"]]*If[p["EdgeLaw"]=="Uniform",r1[[b,i,j]],(r1[[b,i,j]]+r2[[b,i,j]])/2];
    q=If[p["Geometry"]=="Local",If[d==1,q,0],q*adj[[i,j]]/exact[p["DecayBase"]]^(d-If[p["DistanceNormalization"]=="NearestNeighbor",1,0])];
    blocks[[b,i,j]]=q;blocks[[b,j,i]]=q,{b,6},{i,1,k-1},{j,i+1,k}];
  blocks
];
TSGenerateBlocks[p_Association] := Catch[validate[p];generate[p],tsFailure];

TSMixingAngles[u_?MatrixQ] := Module[{a=Abs[u],atan},
  If[Dimensions[u]!={3,3},Return[Failure["Dimensions",<|"MessageTemplate"->"Expected a 3 by 3 active block."|>]]];
  atan[x_,y_] := If[TrueQ[x==0 && y==0],Missing["UndefinedAngle"],ArcTan[x,y]/Degree];
  {atan[a[[1,1]],a[[1,2]]],
   If[TrueQ[a[[1,3]]<=1+10^-20],ArcSin[Min[1,a[[1,3]]]]/Degree,Missing["InvalidActiveNorm"]],
   atan[a[[3,3]],a[[2,3]]]} 
];

solve[blocks_,p_] := Module[
  {k=p["Sites"],n,pr=p["WorkingPrecision"],h,p1,pn,a,b,c,q,rawQ,
   lu,dd,rv,z33,z3n,znn,m,vals,vecs,left,right,signed,order,perm,
   allMasses,lightMasses,lightVectors,rightLight,u,angles,weights,rowWeights,
   residuals,relativeResiduals,degenerate,pairGaps,nonunitary,result,phases,
   symmetryChange,ev,lightPrecision,sourceMachine},
  n=3 k;
  need[Dimensions[blocks]=={6,k,k} && And@@(MatrixQ[#,realQ]& /@ blocks),"Blocks","Supply six real K by K blocks."];
  sourceMachine=!FreeQ[{blocks,Values[p]},_?MachineNumberQ];
  need[And@@(TrueQ[exact[#]==Transpose[exact[#]]]& /@ blocks),"Blocks","Each supplied tower block must be symmetric."];
  h=ArrayFlatten[Table[exact[p["HamiltonianFlavor"][[i,j]]]*exact[blocks[[$blockMap[[i,j]]]]],{i,3},{j,3}]];
  p1=Table[UnitVector[n,1+(i-1) k],{i,3}];pn=Table[UnitVector[n,i k],{i,3}];
  z33=ConstantArray[0,{3,3}];z3n=ConstantArray[0,{3,n}];znn=ConstantArray[0,{n,n}];
  If[p["Nature"]=="Dirac",
    a=exact[p["ActiveCoupling"]]*exact[p["ActiveFlavor"]].p1;
    b=exact[p["DiracRightCoupling"]]*Transpose[pn].Transpose[exact[p["DiracRightFlavor"]]];
    If[p["Method"]=="Rotated",
      {lu,dd,rv}=SingularValueDecomposition[N[h,pr],Tolerance->0];a=N[a,pr].rv;b=Transpose[lu].N[b,pr];h=dd];
    m=N[ArrayFlatten[{{z33,a},{b,h}}],pr];
    {left,dd,right}=SingularValueDecomposition[m,Tolerance->0];
    vals=Diagonal[dd];order=Ordering[vals];allMasses=vals[[order]];
    lightVectors=left[[All,Take[order,3]]];rightLight=right[[All,Take[order,3]]];
    lightMasses=Take[allMasses,3];
    residuals=Table[Max[maxAbs[m.rightLight[[All,j]]-lightMasses[[j]]*lightVectors[[All,j]]],
      maxAbs[Transpose[m].lightVectors[[All,j]]-lightMasses[[j]]*rightLight[[All,j]]]],{j,3}];
    symmetryChange=False,
    a=exact[p["ActiveCoupling"]]*exact[p["ActiveFlavor"]].pn;
    c=exact[p["PsiCoupling"]]*exact[p["PsiFlavor"]].p1;
    rawQ=exact[p["MajoranaScale"]]*exact[p["MajoranaFlavor"]];
    q=If[p["SymmetrizeMajorana"],(rawQ+Transpose[rawQ])/2,rawQ];symmetryChange=!TrueQ[rawQ==q];
    If[p["Method"]=="Rotated",
      {lu,dd,rv}=SingularValueDecomposition[N[h,pr],Tolerance->0];a=N[a,pr].rv;c=N[c,pr].lu;h=dd];
    m=N[ArrayFlatten[{{z33,z33,z3n,a},{z33,q,c,z3n},
      {Transpose[z3n],Transpose[c],znn,h},{Transpose[a],Transpose[z3n],Transpose[h],znn}}],pr];
    need[TrueQ[m==Transpose[m]],"MajoranaSymmetry","Internal Majorana matrix failed symmetry check."];
    {signed,vecs}=Eigensystem[m];
    need[maxAbs[Im[signed]]==0,"Spectrum","Expected real eigenvalues for real symmetric input."];
    vals=Abs[signed];order=Ordering[vals];allMasses=vals[[order]];
    lightMasses=Take[allMasses,3];
    phases=If[#<0,I,1]& /@ signed[[Take[order,3]]];
    lightVectors=Transpose[vecs[[Take[order,3]]]] . DiagonalMatrix[phases];
    residuals=Table[maxAbs[m.lightVectors[[All,j]]-lightMasses[[j]]*Conjugate[lightVectors[[All,j]]]],{j,3}]
  ];
  perm=If[p["MassOrdering"]=="Ascending",{1,2,3},{2,3,1}];
  lightMasses=lightMasses[[perm]];lightVectors=lightVectors[[All,perm]];residuals=residuals[[perm]];
  u=lightVectors[[1;;3,All]];
  relativeResiduals=MapThread[If[TrueQ[#2>0],#1/#2,Missing["ZeroMass"]]&,{residuals,lightMasses}];
  pairGaps=Flatten[Table[If[Max[lightMasses[[i]],lightMasses[[j]]]>0,
    Abs[lightMasses[[i]]-lightMasses[[j]]]/Max[lightMasses[[i]],lightMasses[[j]]],0],{i,1,2},{j,i+1,3}]];
  degenerate=Min[pairGaps]<=p["DegeneracyTolerance"];
  angles=TSMixingAngles[u];weights=Re[Diagonal[ConjugateTranspose[u].u]];
  rowWeights=Re[Diagonal[u.ConjugateTranspose[u]]];
  nonunitary=maxAbs[u.ConjugateTranspose[u]-IdentityMatrix[3]];
  ev=(10^12 lightMasses)^2;
  lightPrecision=Min[Precision /@ lightMasses];
  result=<|"PackageVersion"->$version,"WolframVersion"->$Version,"Parameters"->p,"MatrixDimension"->Length[m],
    "InputIncludedMachineNumbers"->sourceMachine,"MajoranaFlavorWasSymmetrized"->symmetryChange,
    "LightMassesTeV"->lightMasses,"LightMassesEV"->10^12 lightMasses,"AllMassesTeV"->allMasses,
    "PMNS"->u,"AnglesDegrees"->If[degenerate,Missing["DegenerateMasses"],angles],"RawAnglesDegrees"->angles,
    "MassSquaredSplittingsEV2"-><|"Delta21"->ev[[2]]-ev[[1]],"Delta31"->ev[[3]]-ev[[1]],"Delta32"->ev[[3]]-ev[[2]]|>,
    "ActiveColumnWeights"->weights,"ActiveRowWeights"->rowWeights,"MaxNonUnitarity"->nonunitary,
    "FourthMassTeV"->allMasses[[4]],"LightRelativeResiduals"->relativeResiduals,
    "MinimumLightMassPrecision"->lightPrecision,"DegenerateLightMasses"->degenerate,
    "NumericallyResolved"->TrueQ[VectorQ[relativeResiduals,NumericQ] && Max[relativeResiduals]<p["ResidualTolerance"] && lightPrecision>=10]|>;
  If[p["Nature"]=="Majorana",result=Join[result,<|"EffectiveMajoranaBlock"->q|>]];
  If[p["KeepMatrices"],result=Join[result,<|"TowerBlocks"->blocks,"FinalMatrix"->m,"LightVectors"->lightVectors|>]];
  result
];

TSMassFromBlocks[blocks_List,p_Association] := Catch[validate[p];solve[blocks,p],tsFailure];
TSMass[p_Association] := Catch[validate[p];solve[generate[p],p],tsFailure];

Options[TSScan]={"Parallel"->False};
TSScan[p_Association,count_Integer,OptionsPattern[]] := Catch[Module[{jobs,file=$packageFile},
  validate[p];need[count>0,"Count","Sample count must be positive."];
  need[MemberQ[{True,False},OptionValue["Parallel"]],"Parallel","Parallel must be True or False."];
  jobs=Table[Join[p,<|"Seed"->p["Seed"]+i-1|>],{i,count}];
  If[TrueQ[OptionValue["Parallel"]],
    LaunchKernels[];With[{path=file},ParallelEvaluate[Get[path]]];
    ParallelMap[TheorySpaceNeutrinos`TSMass,jobs,DistributedContexts->None],
    TheorySpaceNeutrinos`TSMass /@ jobs]
],tsFailure];

TSPrecisionCheck[p_Association,extraDigits_Integer:20] := Catch[Module[{blocks,a,b,am,bm,ae,me},
  validate[p];need[extraDigits>0,"Precision","extraDigits must be positive."];
  blocks=generate[p];a=solve[blocks,p];b=solve[blocks,Join[p,<|"WorkingPrecision"->p["WorkingPrecision"]+extraDigits|>]];
  am=a["LightMassesTeV"];bm=b["LightMassesTeV"];
  me=If[Min[bm]>0,Max[Abs[am/bm-1]],Missing["ZeroMass"]];
  ae=If[TrueQ[a["DegenerateLightMasses"] || b["DegenerateLightMasses"]],Missing["DegenerateMasses"],maxAbs[a["RawAnglesDegrees"]-b["RawAnglesDegrees"]]];
  <|"LowerPrecision"->p["WorkingPrecision"],"HigherPrecision"->p["WorkingPrecision"]+extraDigits,
    "MaxMassRelativeError"->me,"MaxAngleErrorDegrees"->ae,
    "Pass"->TrueQ[NumericQ[me] && me<10^-8 && NumericQ[ae] && ae<10^-8 && a["NumericallyResolved"] && b["NumericallyResolved"]],
    "LowerResult"->a,"HigherResult"->b|>
],tsFailure];

tableRow[r_Association] := Join[KeyTake[r["Parameters"],{"Seed","Nature","Geometry","Preset","Channel","SiteHalfWidth","HoppingHalfWidth"}],
  <|"LightMassesTeV"->N[r["LightMassesTeV"],10],"AnglesDegrees"->N[r["AnglesDegrees"],10],
    "ActiveWeights"->N[r["ActiveColumnWeights"],8],"MaxNonUnitarity"->N[r["MaxNonUnitarity"],5],
    "NumericallyResolved"->r["NumericallyResolved"],"Degenerate"->r["DegenerateLightMasses"]|>];
tableRow[r_] := <|"Failure"->r|>;
TSTable[results_List] := Dataset[tableRow /@ results];
End[];
EndPackage[];
