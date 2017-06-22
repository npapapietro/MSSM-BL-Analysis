(* ::Package:: *)

BeginPackage["RGEnerate`"]


(* ::Title:: *)
(*The MSSM X Subscript[U(1), B-L] RGEs with U(1) mixings.*)


(* ::Subchapter:: *)
(*Gauge couplings*)


(*General gauge mixing matrix*)
(*Number of U(1)'s*)
Flavors = 6;
Dim = 2;
t0 = N@Log[1000];
t1= Log[1.26*10^16];
(*Use tan\[Beta]=10*)
b = ArcTan[10];
vu = 246*Sin[b];
vd = 246*Cos[b];
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
NN={{Sqrt[3/5],0},{0,Sqrt[2/3]}};

(*U1 Mixing*)
\[Gamma][f_]:= (#[[1,1]]*#[[2,1]]*#[[3,1]]*Take[#,-Dim].Transpose[Take[#,-Dim]])&[f]
gamma = NN.Total[\[Gamma][#]&/@Qs].NN;
G = Table[Symbol["g"<>{ToString[i],ToString[j]} ][t],{i,2},{j,2}];
(*The actual equation is 16 \[Pi]^2 dG/dt=(G)(G^T)(\[Gamma])(G)*)
RHS = Flatten[G.Transpose@G.gamma.G];
LHS = Flatten@G;
G0 = {{0.358839Sqrt[5/3],0},{0,0.1}};
Gi[a_,b_]:=#1==#2/.t->t0&[a,b];
f[x_,y_]:=16 \[Pi]^2  D[#1,t]==#2&[x,y];

{Sg} = NDSolve[Join[MapThread[f,{LHS,RHS}],MapThread[Gi,{LHS,Flatten@G0}]],LHS,{t,t0,t1}];
G1plot =Plot[Evaluate[LHS/. Sg],{t,t0,t1},PlotLabels->{"gyy","gby","gyb","gbb"}];

(*Non Abelian RGEs*)
Gn[a_,b_]:=(16\[Pi]^2*D[#1,t]==#2*(#1)^3)&[a,b];
funcsN = {g2[t],g3[t]};
gammaN = {1,-3};
gN0 = {0.639084,1.06197};
{SgN}=NDSolve[Join[MapThread[Gn,{funcsN,gammaN}],MapThread[Gi,{funcsN,gN0}]],funcsN,{t,t0,t1}];

GNplot = Plot[Evaluate[funcsN/. SgN],{t,t0,t1},PlotLabels->{"g2","g3"},PlotRange->{{t0,t1},{0,2}}];
Show[G1plot,GNplot,PlotRange->{{t0,t1},{0,2}},PlotLabel->"Gauge Couplings"]
(*Remember that gyy is the not the same as g1 or gy from the MSSM, it is in a different basis. That is why it does not show unification in the combined plot*)


(* ::Subchapter:: *)
(*Yukawa couplings*)


(*Yukawa RGEs*)
gammaY[a_]:=(2 Transpose[NN.Take[a,-Dim]].G.Transpose[G].NN.Take[a,-Dim])[[1,1]]

Yu=Table[Symbol["Yu"<>{ToString[i],ToString[j]} ][t],{i,3},{j,3}];
Yd=Table[Symbol["Yd"<>{ToString[i],ToString[j]} ][t],{i,3},{j,3}];
Ye=Table[Symbol["Ye"<>{ToString[i],ToString[j]} ][t],{i,3},{j,3}];
Yui = (2/vu^2){
{2.0653737758081004`*^-6- 7.325449579926939`*^-23 I,7.621443020256718`*^-7- 8.963866581465222`*^-7 I,-0.03263437299165543`+ 0.036974148597956856` I},
{7.62144302025672`*^-7+ 8.963866581465211`*^-7 I,0.3351156380568677` -6.1202969198215815`*^-27 I,-0.5990451797404364`- 5.23233814788266`*^-6 I},
{-0.03263437299165544`- 0.03697414859795678` I,-0.5990451797404364`+ 5.23233814788266`*^-6 I,24136.45377485772` -2.0052964353145634`*^-23 I}};
Ydi =(2/vd^2){
{0.00032130902871975106` -9.05557980138124`*^-24 I,0.0016718956219257289` -0.0007445776642253052` I,0.015687505179799787` -0.017773614156435483` I},
{0.0016718956219257286` +0.0007445776642253052` I,0.016513669322234456` -3.9258029704465054`*^-26 I,0.2879592558345305` +2.512842745645849`*^-6 I},
{0.015687505179799787` +0.017773614156435486` I,0.28795925583453047` -2.512842745645849`*^-6 I,6.878123560897701` -8.631117682651736`*^-24 I}};
Yei = (2/vd)^2{
{2.514728154500983`*^-7- 4.830770091116193`*^-33 I,0,0},
{0,0.010751266003149759` -2.06531025840418`*^-28 I,0},
{0,0,3.0418811625410958` -5.844836374470045`*^-26 I}};
Yui[[3,3]]

g11[t_]:=Evaluate[g11[t]/.Sg];
g12[t_]:=Evaluate[g12[t]/.Sg];
g21[t_]:=Evaluate[g21[t]/.Sg];
g22[t_]:=Evaluate[g22[t]/.Sg];
g2[t_]:=Evaluate[g2[t]/.SgN];
g3[t_]:=Evaluate[g3[t]/.SgN];
\[Beta]u=Yu.(3*ConjugateTranspose[Yu].Yu+ConjugateTranspose[Yd].Yd+3*Tr[ConjugateTranspose[Yu].Yu]-(Total[gammaY[#]&/@{Qq,Qu,QHu}]+3 g2[t]^2+16/3*g3[t]^2)*IdentityMatrix[Flavors/2]);
\[Beta]d=Yd.(3*ConjugateTranspose[Yd].Yd+ConjugateTranspose[Yu].Yu+Tr[3*ConjugateTranspose[Yd].Yd+ConjugateTranspose[Ye].Ye]-(Total[gammaY[#]&/@{Qq,Qd,QHd}]+3 g2[t]^2+16/3 g3[t]^2)*IdentityMatrix[Flavors/2]);
\[Beta]e=Ye.(3*ConjugateTranspose[Ye].Ye+Tr[3*ConjugateTranspose[Yd].Yd+ConjugateTranspose[Ye].Ye]-(Total[gammaY[#]&/@{Ql,Qe,QHd}]+3 g2[t]^2)*IdentityMatrix[Flavors/2]);

Yi[a_,b_]:=#1==Sqrt[#2]/.t->t0&[a,b];
Clear[LHS,RHS]
(* Takes the above Beta functions for Up, Down and Electron Yukawas functions, and turns them all into a list. Then passes them with initial conditions to NDSolve *)
RHS = Join[Flatten@\[Beta]u,Flatten@\[Beta]d,Flatten@\[Beta]e];
LHS = Join[Flatten@Yu,Flatten@Yd,Flatten@Ye];
{Sy} = NDSolve[Join[MapThread[f,{LHS,RHS}],MapThread[Yi,{LHS,Join[Flatten@Yui,Flatten@Ydi,Flatten@Yei]}]],LHS,{t,t0,t1}];
Yup[x_]:=ArrayReshape[Take[Evaluate[LHS/.Sy],9],{3,3}]/.t->x;
Ydo[x_]:=ArrayReshape[Take[Evaluate[LHS/.Sy],{10,18}],{3,3}]/.t->x;
Yel[x_]:=ArrayReshape[Take[Evaluate[LHS/.Sy],-9],{3,3}]/.t->x;
Plot[{Abs[Yup[t][[3,3]]],Abs[Yel[t][[3,3]]]},{t,t0,t1},PlotLabel->"Yukawas"]
Abs[Yup[t1][[3,3]]]


(* ::Subchapter:: *)
(*Gaugino masses*)


(*Gaugino RGEs*)
Mg = Table[Symbol["M"<>{ToString[i],ToString[j]} ][t],{i,2},{j,2}];
Mi[a_,b_]:=#1==#2/.t->t1&[a,b];
\[Beta]M = Mg.Transpose[G].gamma.G +Transpose[G].gamma.G.Mg;

RHSM = Join[Flatten@\[Beta]M,{M2[t],M3[t]}];
LHSM= Join[Flatten@Mg,{M2[t],M3[t]}];
{Sm} = NDSolve[Join[MapThread[f,{LHSM,RHSM}],MapThread[Mi,{LHSM,Table[1000,Length[LHSM]]}]],LHSM,{t,t1,t0}];
Plot[Evaluate[LHSM/.Sm],{t,t1,t0},PlotLabel->"Gaugino Masses",PlotLegends->Table[ToString@i,{i,LHSM}]]



(* ::Subchapter:: *)
(*Soft Masses*)


(* Soft Masses *)
M11[t_]:=Evaluate[M11[t]/.Sm];
M12[t_]:=Evaluate[M12[t]/.Sm];
M21[t_]:=Evaluate[M21[t]/.Sm];
M22[t_]:=Evaluate[M22[t]/.Sm];
M2[t_]:=Evaluate[M2[t]/.Sm];
M3[t_]:=Evaluate[M3[t]/.Sm];

(*Extract Yukawas from NDSolve *)



(*gauge contribution*)
m2g[i_]:=-Transpose[Take[i,-Dim]].G.Mg.ConjugateTranspose[Mg].Transpose[G].Take[i,-Dim];

(*Soft masses functions explicitly*)
Msqe = {Msqe1[t],Msqe2[t],Msqe3[t]};
MsqL = {MsqL1[t],MsqL2[t],MsqL3[t]};
MsqQ = {MsqQ1[t],MsqQ2[t],MsqQ3[t]};
Msqu = {Msqu1[t],Msqu2[t],Msqu3[t]};
Msqd = {Msqd1[t],Msqd2[t],Msqd3[t]};
MsqN = {MsqN1[t],MsqN2[t],MsqN3[t]};
MsqH1 = {{MsqHu[t]}};
MsqH2= {{MsqHd[t]}};
Msqlist = {MsqQ,Msqu,Msqd,MsqL,Msqe,MsqN,MsqH1,MsqH2};

mult[a_,b_]:=#1*#2&[a,b];
traceterm[a_]:=Total[Total[#]&/@MapThread[mult,{Total[Plus@@#]&/@Msqlist,(Transpose[Transpose[G].Take[a,-Dim]].Transpose[G].Take[#,-2]&/@Qs)}]][[1]];
Sum[3*Abs[Yup[t][[i,j]]]^2*(MsqH1[[1]]+MsqQ[[i]]+Msqu[[j]]),{i,3},{j,3}];


(*Define Beta Functions*)
\[Beta]MHu= 1/(8\[Pi]^2)(3*Abs[Yup[t][[3,3]]]^2*(MsqH1[[1]]+MsqQ[[3]]+Msqu[[3]])+traceterm[QHu]+ 8 m2g[QHu]-3 g2[t]^2 M2[t]^2);
\[Beta]MHd = 1/(8\[Pi]^2)(Sum[Abs[Yel[t][[i,j]]]^2*(MsqH2[[1]]+MsqL[[i]]+Msqe[[j]])+3*Abs[Ydo[t][[i,j]]]^2*(MsqH2[[1]]+MsqQ[[i]]+Msqd[[j]]),{i,1,3},{j,1,3}]+traceterm[QHd]+8 m2g[QHd]-3 g2[t]^2 M2[t]^2);
\[Beta]Me[i_]:= 1/(8\[Pi]^2)(Sum[2*Abs[Yel[t][[i,j]]]^2*(MsqH2[[1]]+MsqL[[j]]+Msqe[[i]]),{j,1,3}]+traceterm[Qe]+8 m2g[Qe]);
\[Beta]ML[i_]:=1/(8\[Pi]^2)(Sum[Abs[Yel[t][[i,j]]]^2*(MsqH2[[1]]+MsqL[[i]]+Msqe[[j]]),{j,1,3}]+traceterm[Ql]+8 m2g[Ql]-3 g2[t]^2 M2[t]^2);
\[Beta]Md[i_]:=1/(8\[Pi]^2)(Sum[2*Abs[Ydo[t][[i,j]]]^2*(MsqH2[[1]]+Msqd[[i]]+MsqQ[[j]]),{j,1,3}]+traceterm[Qd]+8 m2g[Qd]-16/3 g3[t]^2 M3[t]^2);
\[Beta]Mu[i_]:=1/(8\[Pi]^2)(Sum[2*Abs[Yup[t][[i,j]]]^2*(MsqH1[[1]]+Msqu[[i]]+MsqQ[[j]]),{j,1,3}]+traceterm[Qu]+8 m2g[Qu]-16/3 g3[t]^2 M3[t]^2);
\[Beta]MQ[i_]:=1/(8\[Pi]^2)(Sum[Abs[Yup[t][[i,j]]]^2*(MsqH1[[1]]+Msqu[[j]]+MsqQ[[i]])+Abs[Ydo[t][[i,j]]]^2*(MsqH2[[1]]+MsqQ[[i]]+Msqd[[j]]),{j,1,3}]+traceterm[Qq]+8 m2g[Qq]-3 g2[t]^2 M2[t]^2-16/3 g3[t]^2 M3[t]^2);
\[Beta]MN[i_]:= 1/(8\[Pi]^2)(traceterm[QN]+8 m2g[QN]); (* Neutrino Yukawas are ignored because of their extreme smallness *)

\[Beta]MHu

LHSMs = Flatten@Msqlist;
RHSMs = Join[Table[\[Beta]MQ[i][[1]],{i,3}],Table[\[Beta]Mu[i][[1]],{i,3}],Table[\[Beta]Md[i][[1]],{i,3}],Table[\[Beta]ML[i][[1]],{i,3}],Table[\[Beta]Me[i][[1]],{i,3}],Table[\[Beta]MN[i][[1]],{i,3}],\[Beta]MHu[[1]],\[Beta]MHd[[1]]];
dterm[x_,y_]:=D[#1,t]==#2[[1]]&[x,y];
initterm[a_,b_]:=#1==#2^2/.t->t1&[a,b];
{Smsq} = NDSolve[Join[MapThread[dterm,{LHSMs,RHSMs}],MapThread[initterm,{LHSMs,Table[1000,Length[LHSMs]]}]],LHSMs,{t,t1,t0}];
Plot[Evaluate[LHSMs/.Smsq],{t,t1,t0},PlotLabel->"Soft Masses",PlotLegends->Table[ToString@i,{i,LHSMs}]]
Plot[Evaluate[Take[LHSMs,-2]/.Smsq],{t,t1,t0},PlotLabel->"Soft Masses",PlotLegends->Table[ToString@i,{i,LHSMs}]]


Exit[]
EndPackage[]
