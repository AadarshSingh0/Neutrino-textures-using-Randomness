(* Independent stored inputs: no Python installation or external paths needed. *)
Get[FileNameJoin[{DirectoryName[$InputFileName],"TheorySpaceNeutrinos.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName],"UnifiedReferenceData.wl"}]];
ClearAll[TSRunChecks,TSGeneratorChecks,TSParallelCheck,tsCheckOne];

tsCheckOne[case_Association] := Module[{p,a,b,ref,me,ae,bme,bae,we,pass},
  p=case["Parameters"];Print["Checking ",case["ID"]];
  a=TSMassFromBlocks[case["Blocks"],Join[p,<|"Method"->"Rotated"|>]];
  b=TSMassFromBlocks[case["Blocks"],Join[p,<|"Method"->"Direct"|>]];
  If[FailureQ[a] || FailureQ[b],Return[<|"Case"->case["ID"],"Check"->"FAIL","Details"->{a,b}|>]];
  ref=case["Reference"];
  me=Max[Abs[a["LightMassesTeV"]/ref["Masses"]-1]];
  ae=Max[Abs[a["RawAnglesDegrees"]-ref["Angles"]]];
  bme=Max[Abs[a["LightMassesTeV"]/b["LightMassesTeV"]-1]];
  bae=Max[Abs[a["RawAnglesDegrees"]-b["RawAnglesDegrees"]]];
  we=Max[Abs[a["ActiveColumnWeights"]-ref["ActiveWeights"]]];
  pass=TrueQ[Max[me,ae,bme,bae,we]<10^-8 && a["NumericallyResolved"] && b["NumericallyResolved"]];
  <|"Case"->case["ID"],"Dimension"->a["MatrixDimension"],
    "MassRelativeErrorVsPython"->N[me,5],"AngleErrorVsPythonDegrees"->N[ae,5],
    "RotatedVsDirectMassRelativeError"->N[bme,5],"RotatedVsDirectAngleErrorDegrees"->N[bae,5],
    "ActiveWeightErrorVsPython"->N[we,5],"Check"->If[pass,"PASS","FAIL"]|>
];

TSRunChecks[which_String:"Smoke"] := Module[{cases},
  cases=Switch[which,
    "Smoke",Select[tsReferenceCases,MemberQ[{
      "Majorana-WeakSite-Local-Hamiltonian","Majorana-StrongSite-Local-Boundary",
      "Majorana-WeakHopping-Local-Hamiltonian","Dirac-WeakSite-Local-Hamiltonian-c3"},#["ID"]]&],
    "All",tsReferenceCases,
    "Prior",Select[tsReferenceCases,TrueQ[#["PreviouslyConfirmed"]]&],
    "StrongHopping",Select[tsReferenceCases,#["Parameters"]["Preset"]=="StrongHopping"&],
    _,Select[tsReferenceCases,#["ID"]==which&]];
  If[cases=={},Return[Failure["Case",<|"MessageTemplate"->"Use Smoke, Prior, All, StrongHopping, or an exact case ID."|>]]];
  tsCheckOne /@ cases
];

TSGeneratorChecks[] := Module[{rows={},record,p,q,b,d,k=8,den,sites,r1,r2,base,expected,adj,off,
  gen,siteScale,law,geom,rn,nu,probe,swap,io,maj,rawQ,asym},
  record[name_,value_] := AppendTo[rows,<|"Test"->name,"Check"->If[TrueQ[value],"PASS","FAIL"]|>];
  p=Join[TSParameters["WeakSite"],<|"Sites"->k,"Seed"->701,"RandomDigits"->12|>];
  b=TSGenerateBlocks[p];record["Same seed reproduces exact blocks",b===TSGenerateBlocks[p]];
  q=Join[p,<|"SiteHalfWidth"->0|>];d=TSGenerateBlocks[q];
  record["Zero disorder removes seed dependence",d===TSGenerateBlocks[Join[q,<|"Seed"->999|>]]];
  q=Join[p,<|"SiteHalfWidth"->2*p["SiteHalfWidth"]|>];d=TSGenerateBlocks[q];
  record["Delta changes only onsite fluctuations",And@@Flatten[Table[
    If[i==j,d[[z,i,j]]-p["SiteMean"]==2*(b[[z,i,j]]-p["SiteMean"]),d[[z,i,j]]==b[[z,i,j]]],{z,6},{i,k},{j,k}]]];
  adj=TSPetersenAdjacency[16];record["Paper Petersen: 16 vertices, 20 edges",Dimensions[adj]=={16,16} && Total[Flatten[adj]]/2==20 && Sort[DeleteDuplicates[Total /@ adj]]=={2,3}];
  record["Petersen adjacency symmetric",adj==Transpose[adj]];
  (* Independently reconstruct generator draws and each expected block entry. *)
  den=10^p["RandomDigits"];
  BlockRandom[SeedRandom[p["Seed"],Method->"MersenneTwister"];
    sites=RandomInteger[{-den,den},{6,k}]/den;
    r1=RandomInteger[{-den,den},{6,k,k}]/den;
    r2=RandomInteger[{-den,den},{6,k,k}]/den];
  Do[q=Join[p,<|"Geometry"->geom,"HoppingMean"->1/7,"HoppingHalfWidth"->2/9,"EdgeLaw"->law|>];
    gen=TSGenerateBlocks[q];adj=If[geom=="Petersen",TSPetersenAdjacency[k],ConstantArray[1,{k,k}]];
    expected=Table[If[i==j,p["SiteMean"]+p["SiteHalfWidth"]*sites[[z,i]],
      base=1/7+(2/9)*If[law=="Uniform",r1[[z,Min[i,j],Max[i,j]]],(r1[[z,Min[i,j],Max[i,j]]]+r2[[z,Min[i,j],Max[i,j]]])/2];
      If[geom=="Local",If[Abs[i-j]==1,base,0],adj[[i,j]]*base/p["DecayBase"]^(Abs[i-j]-1)]],{z,6},{i,k},{j,k}];
    record[geom<>" / "<>law<>" exact generator entries",gen===expected];
    If[geom!="Local",d=TSGenerateBlocks[Join[q,<|"DistanceNormalization"->"Numerator"|>]];
      record[geom<>" / "<>law<>" distance normalization",And@@Flatten[Table[If[i==j,d[[z,i,j]]==gen[[z,i,j]],p["DecayBase"]*d[[z,i,j]]==gen[[z,i,j]]],{z,6},{i,k},{j,k}]]]],
    {geom,{"Local","NonLocal","Petersen"}},{law,{"Uniform","LegacyAverage"}}];
  record["Angle extraction retains upper octant",Abs[TSMixingAngles[{{1,0,0},{0,1/2,Sqrt[3]/2},{0,-Sqrt[3]/2,1/2}}][[3]]-60]<10^-20];
  record["Nonsymmetric Majorana flavor rejected when symmetrization disabled",
    FailureQ[TSGenerateBlocks[Join[TSParameters["WeakSite","Local","Majorana","Boundary"],<|"SymmetrizeMajorana"->False|>]]]];
  record["Unknown parameter rejected",FailureQ[TSGenerateBlocks[Join[p,<|"SiteHalfWidht"->1|>]]]];
  record["Unspecified strong hopping rejected",FailureQ[TSMass[TSParameters["StrongHopping"]]]];
  (* A small fully coupled case verifies ordering and Majorana symmetry. *)
  probe=Join[TSParameters["Custom","Local","Majorana","Boundary"],<|"Sites"->3,"KeepMatrices"->True|>];
  maj=TSMass[probe];
  If[FailureQ[maj],record["Small Majorana construction",False],
    record["Majorana dimension and exact symmetry",Dimensions[maj["FinalMatrix"]]=={24,24} && maj["FinalMatrix"]==Transpose[maj["FinalMatrix"]]];
    record["Takagi light-column residual",maj["NumericallyResolved"]];
    io=TSMass[Join[probe,<|"MassOrdering"->"Inverted"|>]];
    record["Mass order and mixing columns permuted together",!FailureQ[io] && Max[Abs[io["LightMassesTeV"]-maj["LightMassesTeV"][[{2,3,1}]]]]<10^-50 && Max[Flatten[Abs[Abs[io["PMNS"]]-Abs[maj["PMNS"][[All,{2,3,1}]]]]]]<10^-50]
  ];
  rows
];

TSParallelCheck[] := Module[{p,serial,parallel},
  p=Join[TSParameters["WeakSite","Local","Dirac","Hamiltonian"],<|"Sites"->3,"Seed"->819,"WorkingPrecision"->60|>];
  serial=TSScan[p,2];parallel=TSScan[p,2,"Parallel"->True];
  If[!ListQ[parallel] || Length[parallel]!=2 || AnyTrue[Join[serial,parallel],FailureQ],Return[<|"Check"->"FAIL","Serial"->serial,"Parallel"->parallel|>]];
  <|"Test"->"Serial/parallel identical seeded results",
    "Check"->If[And@@MapThread[(Max[Abs[#1["LightMassesTeV"]/#2["LightMassesTeV"]-1]]<10^-20 && Max[Abs[#1["RawAnglesDegrees"]-#2["RawAnglesDegrees"]]]<10^-20)&,{serial,parallel}],"PASS","FAIL"]|>
];
