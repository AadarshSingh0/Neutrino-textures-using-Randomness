(* ::Package:: *)

BeginPackage["NeutrinoModules`"]

makeML::usage="makeML[t, W, K] generates a KxK tridiagonal matrix with random site terms.";
makeMNL::usage="makeMNL[t, W, K, l, c] generates a KxK matrix with exponential decay hopping.";
makeMPET::usage="makeMPET[t, W, K, l, c] generates a KxK matrix based on Petersen graph.";
NeutrinoLocalMass::usage="NeutrinoLocalMass[t, W, K, x, y, z, yuk, v, pr] computes neutrino masses and PMNS matrix for local chain.";
NeutrinoNonLocalMass::usage="NeutrinoNonLocalMass[t, W, K, l, c, x, y, z, yuk, v, pr] computes neutrino masses and PMNS matrix for non-local chain.";
NeutrinoPetMass::usage="NeutrinoPetMass[t, W, K, l, c, x, y, z, yuk, v, pr] computes neutrino masses and PMNS matrix for Petersen graph.";
NeutrinoLocalMajoranaMass::usage="NeutrinoLocalMass[t, W, K, x, y, z, yuk, v, pr] computes neutrino masses and PMNS matrix for local chain.";
NeutrinoNonLocalMajoranaMass::usage="NeutrinoNonLocalMass[t, W, K, l, c, x, y, z, yuk, v, pr] computes neutrino masses and PMNS matrix for non-local chain.";
NeutrinoPetMajoranaMass::usage="NeutrinoPetMass[t, W, K, l, c, x, y, z, yuk, v, pr] computes neutrino masses and PMNS matrix for Petersen graph.";

Begin["`Private`"]
(*Shared random onsite potential generator*)
(*generateG[t_,W_,K_]:=RandomReal[{W-t,W+t},K];*)
generateG[t_,W_,K_]:=RandomReal[{-2W-0t/2,2W+0t/2},K];

(*ml matrix:tridiagonal with random onsite+constant hopping*)
makeML[t_,W_,K_]:=Module[{g,Ox},g=generateG[t,W,K];
Ox=ConstantArray[0,{K,K}];
Do[Ox[[i,i]]=g[[i]],{i,1,K}];
Do[Ox[[i+1,i]]=t;Ox[[i,i+1]]=t,{i,1,K-1}];
Rationalize[Ox,Power[10,-20]]];

(*mnl matrix:exponential decay hopping+random onsite*)
makeMNL[t_,W_,K_,l_,c_]:=Module[{g,Ox},g=generateG[t,W,K];
Ox=Table[l/c^Abs[i-j],{i,1,K},{j,1,K}];
Do[Ox[[i,i]]=g[[i]],{i,1,K}];
Rationalize[Ox,Power[10,-20]]];

(*mpet matrix:modified Petersen graph adjacency+random onsite+exponential decay*)
makeMPET[t_,W_,K_,l_,c_]:=Module[{g,mat,p},g=generateG[t,W,K];
p=PetersenGraph[K/2,K/4];
mat=WeightedAdjacencyMatrix[p];
Do[mat[[i,i]]=g[[i]],{i,1,K}];
Do[Do[If[i=!=j,mat[[i,j]]=l*mat[[i,j]]/c^Abs[i-j]],{j,1,K}],{i,1,K}];
Normal[Rationalize[mat,Power[10,-20]]]];
(*Module to compute neutrino masses and mixing matrices*)NeutrinoLocalMass[t_,W_,K_,ham_,yuk_,v_,pr_:30]:=Module[{ml1,ml2,ml3,ml4,ml5,ml6,TotalML,Lmat,evalue,Rmat,L1,L2,L3,R1,R2,R3,Mma,i,Lmix,mass,Rmix,masses,vr,hamr,yukr},vr=Rationalize[v,Power[10,-pr]];
hamr =Rationalize[ham,Power[10,-pr]]; 
yukr =Rationalize[yuk,Power[10,-pr]]; 
ml1=makeML[t,W,K];ml2=makeML[t,W,K];ml3=makeML[t,W,K];ml4=makeML[t,W,K];ml5=makeML[t,W,K];ml6=makeML[t,W,K];
 TotalML=ArrayFlatten[{{hamr[[1,1]]*ml1,hamr[[1,2]]*ml2,hamr[[1,3]]*ml3},{hamr[[2,1]]*ml2,hamr[[2,2]]*ml4,hamr[[2,3]]*ml5},{hamr[[3,1]]*ml3,hamr[[3,2]]*ml5,hamr[[3,3]]*ml6}}];
(*Perform SVD on TotalML with specified precision*){Lmat,evalue,Rmat}=SingularValueDecomposition[SetPrecision[TotalML,pr]];
(*Extract specific rows*)L1=Lmat[[K]];L2=Lmat[[2 K]];L3=Lmat[[3 K]];R1=Rmat[[1]];R2=Rmat[[1+K]];R3=Rmat[[1+2 K]];
Mma=Table[0,{i,1,3 K+3},{j,1,3 K+3}];
i=1;While[i<=3 K,Mma[[i+3,i+3]]=evalue[[i,i]];i++];
 i=1;While[i<=3 K,Mma[[1,3+i]]=vr*R1[[i]];Mma[[3+i,1]]=vr*L1[[i]]*yukr[[1,1]]+vr*L2[[i]]*yukr[[1,2]]+vr*L3[[i]]*yukr[[1,3]];Mma[[2,3+i]]=vr*R2[[i]];Mma[[3+i,2]]=vr*L1[[i]]*yukr[[2,1]]+vr*L2[[i]]*yukr[[2,2]]+vr*L3[[i]]*yukr[[2,3]];Mma[[3,3+i]]=vr*R3[[i]];Mma[[3+i,3]]=vr*L1[[i]]*yukr[[3,1]]+vr*L2[[i]]*yukr[[3,2]]+vr*L3[[i]]*yukr[[3,3]];i++];
{Lmix,mass,Rmix}=SingularValueDecomposition[SetPrecision[Mma,pr]];
masses=Diagonal[mass];PMNS = Table[0,{i,1,3},{j,1,3}];i=1;
While[i<4,{j=1;
While[j<4,{PMNS[[j,-i]]=Lmix[[j,-4+i]]};j++]};
i++] ;
{masses,PMNS}]
(*Module to compute neutrino masses and mixing matrices*)NeutrinoNonLocalMass[t_,W_,K_,l_,c_,ham_,yuk_,v_,pr_:30]:=Module[{ml1,ml2,ml3,ml4,ml5,ml6,TotalML,Lmat,evalue,Rmat,L1,L2,L3,R1,R2,R3,Mma,i,Lmix,mass,Rmix,masses,vr,yukr,hamr},vr=Rationalize[v,Power[10,-pr]];hamr =Rationalize[ham,Power[10,-pr]];yukr =Rationalize[yuk,Power[10,-pr]]; ml1=makeMNL[t,W,K,l,c];ml2=makeMNL[t,W,K,l,c];ml3=makeMNL[t,W,K,l,c];ml4=makeMNL[t,W,K,l,c];ml5=makeMNL[t,W,K,l,c];ml6=makeMNL[t,W,K,l,c];TotalML=ArrayFlatten[{{hamr[[1,1]]*ml1,hamr[[1,2]]*ml2,hamr[[1,3]]*ml3},{hamr[[2,1]]*ml2,hamr[[2,2]]*ml4,hamr[[2,3]]*ml5},{hamr[[3,1]]*ml3,hamr[[3,2]]*ml5,hamr[[3,3]]*ml6}}];
(*Perform SVD on TotalML with specified precision*){Lmat,evalue,Rmat}=SingularValueDecomposition[SetPrecision[TotalML,pr]];
(*Extract specific rows*)L1=Lmat[[K]];L2=Lmat[[2 K]];L3=Lmat[[3 K]];R1=Rmat[[1]];R2=Rmat[[1+K]];R3=Rmat[[1+2 K]];
Mma=Table[0,{i,1,3 K+3},{j,1,3 K+3}];
i=1;While[i<=3 K,Mma[[i+3,i+3]]=evalue[[i,i]];i++];
 i=1;While[i<=3 K,Mma[[1,3+i]]=vr*R1[[i]];Mma[[3+i,1]]=vr*L1[[i]]*yukr[[1,1]]+vr*L2[[i]]*yukr[[1,2]]+vr*L3[[i]]*yukr[[1,3]];Mma[[2,3+i]]=vr*R2[[i]];Mma[[3+i,2]]=vr*L1[[i]]*yukr[[2,1]]+vr*L2[[i]]*yukr[[2,2]]+vr*L3[[i]]*yukr[[2,3]];Mma[[3,3+i]]=vr*R3[[i]];Mma[[3+i,3]]=vr*L1[[i]]*yukr[[3,1]]+vr*L2[[i]]*yukr[[3,2]]+vr*L3[[i]]*yukr[[3,3]];i++];
{Lmix,mass,Rmix}=SingularValueDecomposition[SetPrecision[Mma,pr]];
masses=Diagonal[mass];PMNS = Table[0,{i,1,3},{j,1,3}];i=1;
While[i<4,{j=1;
While[j<4,{PMNS[[j,-i]]=Lmix[[j,-4+i]]};j++]};
i++] ;
{masses,PMNS}]
(*Module to compute neutrino masses and mixing matrices*)NeutrinoPetMass[t_,W_,K_,l_,c_,ham_,yuk_,v_,pr_:30]:=Module[{ml1,ml2,ml3,ml4,ml5,ml6,TotalML,Lmat,evalue,Rmat,L1,L2,L3,R1,R2,R3,Mma,i,Lmix,mass,Rmix,masses,vr,yukr,hamr},vr=Rationalize[v,Power[10,-pr]];yukr =Rationalize[yuk,Power[10,-pr]];hamr =Rationalize[ham,Power[10,-pr]];ml1=makeMPET[t,W,K,l,c];ml2=makeMPET[t,W,K,l,c];ml3=makeMPET[t,W,K,l,c];ml4=makeMPET[t,W,K,l,c];ml5=makeMPET[t,W,K,l,c];ml6=makeMPET[t,W,K,l,c];TotalML=ArrayFlatten[{{hamr[[1,1]]*ml1,hamr[[1,2]]*ml2,hamr[[1,3]]*ml3},{hamr[[2,1]]*ml2,hamr[[2,2]]*ml4,hamr[[2,3]]*ml5},{hamr[[3,1]]*ml3,hamr[[3,2]]*ml5,hamr[[3,3]]*ml6}}];(*Perform SVD on TotalML with specified precision*){Lmat,evalue,Rmat}=SingularValueDecomposition[SetPrecision[TotalML,pr]];
(*Extract specific rows*)L1=Lmat[[K]];L2=Lmat[[2 K]];L3=Lmat[[3 K]];R1=Rmat[[1]];R2=Rmat[[1+K]];R3=Rmat[[1+2 K]];
Mma=Table[0,{i,1,3 K+3},{j,1,3 K+3}];
i=1;While[i<=3 K,Mma[[i+3,i+3]]=evalue[[i,i]];i++];
 i=1;While[i<=3 K,Mma[[1,3+i]]=vr*R1[[i]];Mma[[3+i,1]]=vr*L1[[i]]*yukr[[1,1]]+vr*L2[[i]]*yukr[[1,2]]+vr*L3[[i]]*yukr[[1,3]];Mma[[2,3+i]]=vr*R2[[i]];Mma[[3+i,2]]=vr*L1[[i]]*yukr[[2,1]]+vr*L2[[i]]*yukr[[2,2]]+vr*L3[[i]]*yukr[[2,3]];Mma[[3,3+i]]=vr*R3[[i]];Mma[[3+i,3]]=vr*L1[[i]]*yukr[[3,1]]+vr*L2[[i]]*yukr[[3,2]]+vr*L3[[i]]*yukr[[3,3]];i++];
{Lmix,mass,Rmix}=SingularValueDecomposition[SetPrecision[Mma,pr]];
masses=Diagonal[mass];PMNS = Table[0,{i,1,3},{j,1,3}];i=1;
While[i<4,{j=1;
While[j<4,{PMNS[[j,-i]]=Lmix[[j,-4+i]]};j++]};
i++] ;
{masses,PMNS}]
(*Module to compute neutrino masses and mixing matrices*)NeutrinoLocalMajoranaMass[t_,W_,K_,ham_,yuk_,v_,pr_:30]:=Module[{ml1,ml2,ml3,ml4,ml5,ml6,TotalML,Lmat,evalue,Rmat,L1,L2,L3,R1,R2,R3,Mma,i,i1,j1,Lmix,Ox,Mm2,a1,a2,a3,ac,Mm,ar1,ar2,ar3,Mn,mass,Rmix,masses,PMNS,Mm1,vr,yukr,hamr,Wr},vr=Rationalize[v,Power[10,-pr]];Wr=Rationalize[W,Power[10,-pr]];hamr =Rationalize[ham,Power[10,-pr]];yukr =Rationalize[yuk,Power[10,-pr]]; ml1=makeML[t,W,K];ml2=makeML[t,W,K];ml3=makeML[t,W,K];ml4=makeML[t,W,K];ml5=makeML[t,W,K];ml6=makeML[t,W,K];
TotalML=ArrayFlatten[{{hamr[[1,1]]*ml1,hamr[[1,2]]*ml2,hamr[[1,3]]*ml3},{hamr[[2,1]]*ml2,hamr[[2,2]]*ml4,hamr[[2,3]]*ml5},{hamr[[3,1]]*ml3,hamr[[3,2]]*ml5,hamr[[3,3]]*ml6}}];Ox = Table[0,{i,1,3K},{j,1,3K}]; a1 = Table[0,{i,1,2K*3}];
a2 = Table[0,{i,1,2K*3}];
a3 = Table[0,{i,1,2K*3}];
ac = {a1,a2,a3};i=1;
While[i<4,{ac[[i,1+K*(i-1)]]=t};i++];Mm = {{Ox,TotalML},{Transpose[TotalML],Ox}};
Mm1 =ArrayFlatten[Mm];Mm2 = Transpose[Join[Mm1,{ac[[1]]},{ac[[2]]},{ac[[3]]}]];
ar1 = Append[ac[[1]],yukr[[1,1]]*Wr];AppendTo[ar1,yukr[[1,2]]*Wr];AppendTo[ar1,yukr[[1,3]]*Wr];
ar2 = Append[ac[[2]],yukr[[2,1]]*Wr];AppendTo[ar2,yukr[[2,2]]*Wr];AppendTo[ar2,yukr[[2,3]]*Wr];
ar3 = Append[ac[[3]],yukr[[3,1]]*Wr];AppendTo[ar3,yukr[[3,2]]*Wr];AppendTo[ar3,yukr[[3,3]]*Wr];Mn = Join[Mm2,{ar1},{ar2},{ar3}] ;
(*Perform SVD on TotalML with specified precision*){Lmat,evalue,Rmat}=SingularValueDecomposition[SetPrecision[Mn,pr]];
(*Extract specific rows*)L1=Lmat[[K]];L2=Lmat[[2 K]];L3=Lmat[[3 K]];R1=Rmat[[1+3*K]];R2=Rmat[[1+K+3*K]];R3=Rmat[[1+2 K+3*K]];
Mma=Table[0,{i,1,6 K+6},{j,1,6 K+6}];
i=1;While[i<2K*3+4,Mma[[i+3,i+3]]=evalue[[i,i]];i++];
 i=1;While[i<2K*3+4,Mma[[1,3+i]]=vr*R1[[i]];Mma[[3+i,1]]=vr*L1[[i]];Mma[[2,3+i]]=vr*R2[[i]];Mma[[3+i,2]]=vr*L2[[i]];Mma[[3,3+i]]=vr*R3[[i]];Mma[[3+i,3]]=vr*L3[[i]];i++];
{Lmix,mass,Rmix}=SingularValueDecomposition[SetPrecision[Mma,pr]];
masses=Diagonal[mass];PMNS = Table[0,{i,1,3},{j,1,3}];i1=1;
While[i1<4,{j1=1;
While[j1<4,{PMNS[[j1,-i1]]=Lmix[[j1,-4+i1]]};j1++]};
i1++] ;
{masses,PMNS}]
(*Module to compute neutrino masses and mixing matrices*)NeutrinoNonLocalMajoranaMass[t_,W_,K_,l_,c_,ham_,yuk_,v_,pr_:30]:=Module[{ml1,ml2,ml3,ml4,ml5,ml6,TotalML,Lmat,evalue,Rmat,L1,L2,L3,R1,R2,R3,Mma,i,i1,j1,Lmix,Ox,Mm2,a1,a2,a3,ac,Mm,ar1,ar2,ar3,Mn,mass,Rmix,masses,PMNS,Mm1,vr,yukr,hamr,Wr},vr=Rationalize[v,Power[10,-pr]];Wr=Rationalize[W,Power[10,-pr]];hamr =Rationalize[ham,Power[10,-pr]];yukr =Rationalize[yuk,Power[10,-pr]];ml1=makeMNL[t,W,K,l,c];ml2=makeMNL[t,W,K,l,c];ml3=makeMNL[t,W,K,l,c];ml4=makeMNL[t,W,K,l,c];ml5=makeMNL[t,W,K,l,c];ml6=makeMNL[t,W,K,l,c];
TotalML=ArrayFlatten[{{hamr[[1,1]]*ml1,hamr[[1,2]]*ml2,hamr[[1,3]]*ml3},{hamr[[2,1]]*ml2,hamr[[2,2]]*ml4,hamr[[2,3]]*ml5},{hamr[[3,1]]*ml3,hamr[[3,2]]*ml5,hamr[[3,3]]*ml6}}];
Ox = Table[0,{i,1,3K},{j,1,3K}]; a1 = Table[0,{i,1,2K*3}];
a2 = Table[0,{i,1,2K*3}];
a3 = Table[0,{i,1,2K*3}];
ac = {a1,a2,a3};i=1;
While[i<4,{ac[[i,1+K*(i-1)]]=t};i++];Mm = {{Ox,TotalML},{Transpose[TotalML],Ox}};
Mm1 =ArrayFlatten[Mm];Mm2 = Transpose[Join[Mm1,{ac[[1]]},{ac[[2]]},{ac[[3]]}]];
ar1 = Append[ac[[1]],yukr[[1,1]]*Wr];AppendTo[ar1,yukr[[1,2]]*Wr];AppendTo[ar1,yukr[[1,3]]*Wr];
ar2 = Append[ac[[2]],yukr[[2,1]]*Wr];AppendTo[ar2,yukr[[2,2]]*Wr];AppendTo[ar2,yukr[[2,3]]*Wr];
ar3 = Append[ac[[3]],yukr[[3,1]]*Wr];AppendTo[ar3,yukr[[3,2]]*Wr];AppendTo[ar3,yukr[[3,3]]*Wr];Mn = Join[Mm2,{ar1},{ar2},{ar3}] ;
(*Perform SVD on TotalML with specified precision*){Lmat,evalue,Rmat}=SingularValueDecomposition[SetPrecision[Mn,pr]];
(*Extract specific rows*)L1=Lmat[[K]];L2=Lmat[[2 K]];L3=Lmat[[3 K]];R1=Rmat[[1+3*K]];R2=Rmat[[1+K+3*K]];R3=Rmat[[1+2 K+3*K]];
Mma=Table[0,{i,1,6 K+6},{j,1,6 K+6}];
i=1;While[i<2K*3+4,Mma[[i+3,i+3]]=evalue[[i,i]];i++];
 i=1;While[i<2K*3+4,Mma[[1,3+i]]=vr*R1[[i]];Mma[[3+i,1]]=vr*L1[[i]];Mma[[2,3+i]]=vr*R2[[i]];Mma[[3+i,2]]=vr*L2[[i]];Mma[[3,3+i]]=vr*R3[[i]];Mma[[3+i,3]]=vr*L3[[i]];i++];
{Lmix,mass,Rmix}=SingularValueDecomposition[SetPrecision[Mma,pr]];
masses=Diagonal[mass];PMNS = Table[0,{i,1,3},{j,1,3}];i1=1;
While[i1<4,{j1=1;
While[j1<4,{PMNS[[j1,-i1]]=Lmix[[j1,-4+i1]]};j1++]};
i1++] ;
{masses,PMNS}]
(*Module to compute neutrino masses and mixing matrices*)NeutrinoPetMajoranaMass[t_,W_,K_,l_,c_,ham_,yuk_,v_,pr_:30]:=Module[{ml1,ml2,ml3,ml4,ml5,ml6,TotalML,Lmat,evalue,Rmat,L1,L2,L3,R1,R2,R3,Mma,i,i1,j1,Lmix,Ox,Mm2,a1,a2,a3,ac,Mm,ar1,ar2,ar3,Mn,mass,Rmix,masses,PMNS,Mm1,vr,yukr,hamr,Wr},vr=Rationalize[v,Power[10,-pr]];Wr=Rationalize[W,Power[10,-pr]];hamr =Rationalize[ham,Power[10,-pr]];yukr =Rationalize[yuk,Power[10,-pr]];ml1=makeMPET[t,W,K,l,c];ml2=makeMPET[t,W,K,l,c];ml3=makeMPET[t,W,K,l,c];ml4=makeMPET[t,W,K,l,c];ml5=makeMPET[t,W,K,l,c];ml6=makeMPET[t,W,K,l,c];
TotalML=ArrayFlatten[{{hamr[[1,1]]*ml1,hamr[[1,2]]*ml2,hamr[[1,3]]*ml3},{hamr[[2,1]]*ml2,hamr[[2,2]]*ml4,hamr[[2,3]]*ml5},{hamr[[3,1]]*ml3,hamr[[3,2]]*ml5,hamr[[3,3]]*ml6}}];
Ox = Table[0,{i,1,3K},{j,1,3K}]; a1 = Table[0,{i,1,2K*3}];
a2 = Table[0,{i,1,2K*3}];
a3 = Table[0,{i,1,2K*3}];
ac = {a1,a2,a3};i=1;
While[i<4,{ac[[i,1+K*(i-1)]]=t};i++];Mm = {{Ox,TotalML},{Transpose[TotalML],Ox}};
Mm1 =ArrayFlatten[Mm];Mm2 = Transpose[Join[Mm1,{ac[[1]]},{ac[[2]]},{ac[[3]]}]];
ar1 = Append[ac[[1]],yukr[[1,1]]*Wr];AppendTo[ar1,yukr[[1,2]]*Wr];AppendTo[ar1,yukr[[1,3]]*Wr];
ar2 = Append[ac[[2]],yukr[[2,1]]*Wr];AppendTo[ar2,yukr[[2,2]]*Wr];AppendTo[ar2,yukr[[2,3]]*Wr];
ar3 = Append[ac[[3]],yukr[[3,1]]*Wr];AppendTo[ar3,yukr[[3,2]]*Wr];AppendTo[ar3,yukr[[3,3]]*Wr];Mn = Join[Mm2,{ar1},{ar2},{ar3}] ;
(*Perform SVD on TotalML with specified precision*){Lmat,evalue,Rmat}=SingularValueDecomposition[SetPrecision[Mn,pr]];
(*Extract specific rows*)L1=Lmat[[K]];L2=Lmat[[2 K]];L3=Lmat[[3 K]];R1=Rmat[[1+3*K]];R2=Rmat[[1+K+3*K]];R3=Rmat[[1+2 K+3*K]];
Mma=Table[0,{i,1,6 K+6},{j,1,6 K+6}];
i=1;While[i<2K*3+4,Mma[[i+3,i+3]]=evalue[[i,i]];i++];
 i=1;While[i<2K*3+4,Mma[[1,3+i]]=vr*R1[[i]];Mma[[3+i,1]]=vr*L1[[i]];Mma[[2,3+i]]=vr*R2[[i]];Mma[[3+i,2]]=vr*L2[[i]];Mma[[3,3+i]]=vr*R3[[i]];Mma[[3+i,3]]=vr*L3[[i]];i++];
{Lmix,mass,Rmix}=SingularValueDecomposition[SetPrecision[Mma,pr]];
masses=Diagonal[mass];PMNS = Table[0,{i,1,3},{j,1,3}];i1=1;
While[i1<4,{j1=1;
While[j1<4,{PMNS[[j1,-i1]]=Lmix[[j1,-4+i1]]};j1++]};
i1++] ;
{masses,PMNS}]
End[]
EndPackage[]

