(* ::Package:: *)

BeginPackage["RGEnerate`"]
(*General gauge mixing matrix*)
(*Number of U(1)'s*)
Dim = 2;
t0 = 1000;
t1= 1.26*10^16;


(*Define a vector of particle's charges under SU(N) and U(1)'s--\[Rule] (generation,SU(3),SU(2),Y,B-L)*)
Qq={{3},{3},{2},{1/6},{1/3}};
Qu={{3},{3},{1},{-2/3},{-1/3}};
Qd={{3},{3},{1},{1/3},{-1/3}};
Ql={{3},{1},{2},{-1/2},{-1}};
Qe={{3},{1},{1},{1},{1}};
QN={{3},{1},{1},{0},{1}};
QHu={{1},{1},{2},{1/2},{0}};
QHd={{1},{1},{2},{-1/2},{0}};

Qs = {Qq,Qu,Qd,Ql,Qe,QN,QHu,QHd};
(*Normalization with GUTS ? *)
NN={{1,0},{0,1}};

(*U1 Mixing*)
\[Gamma][f_]:= (#[[1,1]]*#[[2,1]]*#[[3,1]]*Take[#,-Dim].Transpose[Take[#,-Dim]])&[f]
gamma = Total[\[Gamma][#]&/@Qs];
funcs= {{gyy[t],gby[t]},{gyb[t],gbb[t]}};
(*The actual equation is 16 \[Pi]^2 dG/dt=(G)(G^T)(\[Gamma])(G)*)
RHS = Flatten[funcs.Transpose@funcs.gamma.funcs];
LHS = Flatten@funcs;
G0 = {{0.465625,0},{0,0.1}};
Gi[a_,b_]:=#1==#2/.t->t0&[a,b];
f[x_,y_]:=16 \[Pi]^2 t D[#1,t]==#2&[x,y];
Sg = NDSolve[Join[MapThread[f,{LHS,RHS}],MapThread[Gi,{LHS,Flatten@G0}]],funcs,{t,t0,t1}];
Plot[Evaluate[funcs/. Sg],{t,t0,t1},PlotLabels->{"gyy","gby","gyb","gbb"}]


(*SUN RGEs*)
Gn[a_,b_]:=(16\[Pi]^2*D[#1,t]t==#2*(#1)^3)&[a,b];
GNi[a_,b_]:=#1==#2/.t->t0&[a,b];
funcsN = {g2[t],g3[t]};
gammaN = {1,-7};
gN0 = {0.636093,1.06197};
SgN=NDSolve[Join[MapThread[Gn,{funcsN,gammaN}],MapThread[GNi,{funcsN,gN0}]]
,funcsN,{t,t0,t1}]
Plot[Evaluate[funcsN/. SgN],{t,t0,t1},PlotLabels->{"g2","g3"}]


EndPackage[]



