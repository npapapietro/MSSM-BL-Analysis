(* ::Package:: *)

BeginPackage["RGEnerate`"]

(*General gauge mixing matrix*)
(*Number of U(1)'s*)
Dim = 2;
t0 = Log[1000];
t1= Log[1.26*10^16];
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
funcs= {{gyy[t],gby[t]},{gyb[t],gbb[t]}};
(*The actual equation is 16 \[Pi]^2 dG/dt=(G)(G^T)(\[Gamma])(G)*)
RHS = Flatten[funcs.Transpose@funcs.gamma.funcs];
LHS = Flatten@funcs;
G0 = {{0.358839,0},{0,0.1}};
Gi[a_,b_]:=#1==#2/.t->t0&[a,b];
f[x_,y_]:=16 \[Pi]^2  D[#1,t]==#2&[x,y];
{Sg} = NDSolve[Join[MapThread[f,{LHS,RHS}],MapThread[Gi,{LHS,Flatten@G0}]],LHS,{t,t0,t1}];
G1plot =Plot[Evaluate[LHS/. Sg],{t,t0,t1},PlotLabels->{"gyy","gby","gyb","gbb"}];

G = {{gyy[t],gby[t]},{gyb[t],gbb[t]}};

(*SUN RGEs*)
Gn[a_,b_]:=(16\[Pi]^2*D[#1,t]==#2*(#1)^3)&[a,b];

funcsN = {g2[t],g3[t]};
gammaN = {1,-3};
gN0 = {0.639084,1.06197};
{SgN}=NDSolve[Join[MapThread[Gn,{funcsN,gammaN}],MapThread[Gi,{funcsN,gN0}]]
,funcsN,{t,t0,t1}];
GNplot = Plot[Evaluate[funcsN/. SgN],{t,t0,t1},PlotLabels->{"g2","g3"},PlotRange->{{t0,t1},{0,2}}];
Show[G1plot,GNplot,PlotRange->{{t0,t1},{0,2}}]
(*Remember that gyy is the not the same as g1 or gy from the MSSM, it is in a different basis. That is why it does not show unification in the combined plot*)


(*Yukawa RGEs*)
gammaY[a_]:=(2 Transpose[NN.Take[a,-Dim]].G.Transpose[G].NN.Take[a,-Dim])[[1,1]]*IdentityMatrix[3];
(*Use tan\[Beta]=10*)
b = ArcTan[10];
vu = 246*Sin[b];
vd = 246*Cos[b];
Yu={
{Yu11[t],Yu12[t],Yu13[t]},
{Yu21[t],Yu22[t],Yu23[t]},
{Yu31[t],Yu32[t],Yu33[t]}};
Yd={
{Yd11[t],Yd12[t],Yd13[t]},
{Yd21[t],Yd22[t],Yd23[t]},
{Yd31[t],Yd32[t],Yd33[t]}};
Ye={
{Ye11[t],Ye12[t],Ye13[t]},
{Ye21[t],Ye22[t],Ye23[t]},
{Ye31[t],Ye32[t],Ye33[t]}};
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

\[Beta]u=Yu.(3*Conjugate@Transpose@Yu.Yu+Conjugate@Transpose@Yd.Yd+3*Tr[Conjugate@Transpose@Yu.Yu]-Total[gammaY[#]&/@{Qq,Qu,QHu}]-(3 g2[t]^2+16/3 g3[t])^2*IdentityMatrix[3]);
\[Beta]d=Yd.(3*Conjugate@Transpose@Yd.Yd+Conjugate@Transpose@Yu.Yu+Tr[3*Conjugate@Transpose@Yd.Yd+Conjugate@Transpose@Ye.Ye]-Total[gammaY[#]&/@{Qq,Qd,QHd}]-(3 g2[t]^2+16/3 g3[t])^2*IdentityMatrix[3]);
\[Beta]e=Ye.(3*Conjugate@Transpose@Ye.Ye+Tr[3*Conjugate@Transpose@Yd.Yd+Conjugate@Transpose@Ye.Ye]-Total[gammaY[#]&/@{Ql,Qe,QHd}]-(3 g2[t]^2)^2*IdentityMatrix[3]);

Yi[a_,b_]:=#1==Sqrt[#2]/.t->t0&[a,b];
Clear[RHS,LHS];
RHS = Join[Flatten@\[Beta]u,Flatten@\[Beta]d,Flatten@\[Beta]e];
LHS = Join[Flatten@Yu,Flatten@Yd,Flatten@Ye];
gyy[t_]:=Evaluate[gyy[t]/.Sg];
gby[t_]:=Evaluate[gby[t]/.Sg];
gyb[t_]:=Evaluate[gyb[t]/.Sg];
gbb[t_]:=Evaluate[gbb[t]/.Sg];
g2[t_]:=Evaluate[g2[t]/.SgN];
g3[t_]:=Evaluate[g3[t]/.SgN];
{Sy} = NDSolve[Join[MapThread[f,{LHS,RHS}],MapThread[Yi,{LHS,Join[Flatten@Yui,Flatten@Ydi,Flatten@Yei]}]],LHS,{t,t0,t1}];
Plot[Evaluate[LHS/.Sy],{t,t0,t1}]


Exit[]
EndPackage[]



