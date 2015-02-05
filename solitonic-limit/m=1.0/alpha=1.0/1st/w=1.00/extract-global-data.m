(* ::Package:: *)

(*24 Sept 14 *)



(*<<NumericalMath`ListIntegrate` *)


<<NumericalMath`ListIntegrate`
Off[SetDelayed::"write"]


 (*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	              PARTE TOTAL MASS  CALCULATA INTEGRAL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


Directory[]


Off[General::spell1]
Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];

(* 1  2   3   4  5  6   7*)
(* nr,w,alpha,c1,c2,c3,Lambda*)
conf=ReadList["res.txt",{Number,Number ,Number ,Number ,Number,Number  }];

 
nr=conf[[1]][[1]];
wf= conf[[1]][[2]];
alfa= conf[[1]][[3]];
c1= conf[[1]][[4]];
c2= conf[[1]][[5]];
c3= conf[[1]][[6]];  

Print["winding number n = ",nr];
Print["frequency w = ",wf];
Print["alfa = ",alfa];
Print["c1(Phi^6) = ",c1];
Print["c2(Phi^4) = ",c2];
Print["c3(Phi^2) = ",c3]; 


gr=ReadList["gridx.dat",{Number}];


lgr=Length[gr];
nx=lgr;

Print["nx = ",nx];

listar=Table[gr[[k]][[1]],{k,1,lgr}] ;
listalogr=Table[Log[10,gr[[k]][[1]]],{k,1,lgr}];

 unghi0=ReadList["gridy.dat",{Number}];
ny=Length[unghi0];
Print["ny = ",ny];
unghi=Table[unghi0[[k]][[1]],{k,1,ny}]; 

(*unghi=Table[(k-1)*Pi/2/(ny-1),{k,1,ny}];*)

ntot=nx*ny;

a=ReadList["functf.dat",{Number,Number,Number,Number,Number  }];


lung1=Length[a];


(*Datele sunt salvate direct cu indexarea globala*)
f1=Table[a[[k]][[1]],{k,1,lung1}];
f2=Table[a[[k]][[2]],{k,1,lung1}];
f0=Table[a[[k]][[3]],{k,1,lung1}]; 
Z=Table[a[[k]][[4]],{k,1,lung1}]; 
w=Table[a[[k]][[5]],{k,1,lung1}]; 
  
(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[

f1u[k]=Table[f1[[i]],{i,(k-1)*nx+1,k*nx}];
f2u[k]=Table[f2[[i]],{i,(k-1)*nx+1,k*nx}];
f0u[k]=Table[f0[[i]],{i,(k-1)*nx+1,k*nx}]; 
wu[k]=Table[w[[i]],{i,(k-1)*nx+1,k*nx}];
Zu[k]=Table[Z[[i]],{i,(k-1)*nx+1,k*nx}]; 
,{k,1,ny}];
 
 
as1=2;
as2=IntegerPart[ny/2];
as3=ny-1;

sa1=3;
sa2=IntegerPart[nx/2];
sa3=nx-1;


Print["rmax = ",gr[[nx]][[1]]];


minf0=Min[f0] 
f0H=f0u[1][[1]]
f1H=f1u[1][[1]]
maxW=Max[Abs[w]] 
maxZ=Max[Abs[Z]] 



n1=nx-23;

t1=Table[ {listar[[i]] ,E^(2*f0u[as1][[i]])},{i,2,n1}];
t2=Table[{ listar[[i]]  ,E^(2*f0u[as2][[i]])},{i,2,n1}];
t3=Table[{ listar[[i]] ,E^(2*f0u[as3][[i]])},{i,2,n1}];

g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity,  Frame->True ,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ];

g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True ];
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]];

as=Table[{  f0u[as3][[i]]},{i,2,n1}] ;
Min[f0]


t1=Table[{i,f1u[as1][[i]]},{i,1,lgr}];
t2=Table[{i,f1u[as2][[i]]},{i,1,lgr}];
t3=Table[{i,f1u[as3][[i]]},{i,1,lgr}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" FUNCTION f1"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]

t1=Table[ {listar[[i]] ,f1u[as1][[i]]},{i,2,lgr-1}];
t2=Table[{ listar[[i]]  ,f1u[as2][[i]]},{i,2,lgr-1}];
t3=Table[{ listar[[i]] ,f1u[as3][[i]]},{i,2,lgr-1}];

g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity,  Frame->True ,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ];

g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True ];
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]



t1=Table[{i,f2u[as1][[i]]},{i,1,lgr}];
t2=Table[{i,f2u[as2][[i]]},{i,1,lgr}];
t3=Table[{i,f2u[as3][[i]]},{i,1,lgr}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" FUNCTION f2"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]


t1=Table[ {listar[[i]] ,f2u[as1][[i]]},{i,2,lgr-14}];
t2=Table[{ listar[[i]]  ,f2u[as2][[i]]},{i,2,lgr-14}];
t3=Table[{ listar[[i]] ,f2u[as3][[i]]},{i,2,lgr-14}];

g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity,  Frame->True ,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ];

g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True ];
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]


t1=Table[{i,wu[as1][[i]]},{i,1,lgr}];
t2=Table[{i,wu[as2][[i]]},{i,1,lgr}];
t3=Table[{i,wu[as3][[i]]},{i,1,lgr}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" FUNCTION w"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]


t1=Table[ {listar[[i]] ,wu[as1][[i]]},{i,2,lgr-1 }];
t2=Table[{ listar[[i]]  ,wu[as2][[i]]},{i,2,lgr-1 }];
t3=Table[{ listar[[i]] ,wu[as3][[i]]},{i,2,lgr-1 }];

g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity,  Frame->True ,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ];

g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True ];
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]


t1=Table[{i,Zu[as1][[i]] },{i,1,lgr}];
t2=Table[{i,Zu[as2][[i]] },{i,1,lgr}];
t3=Table[{i,Zu[as3][[i]] },{i,1,lgr}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" FUNCTION Z"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]

t1=Table[ {listar[[i]] ,Zu[as1][[i]]},{i,2,lgr-2}];
t2=Table[{ listar[[i]]  ,Zu[as2][[i]]},{i,2,lgr-2}];
t3=Table[{ listar[[i]] ,Zu[as3][[i]]},{i,2,lgr-2}];

g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity,  Frame->True ,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True ];
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]

 







(* now I verify the 1/r decay of F1,F2,F0 *)


ni=100;
t1=Table[{listar[[i]],-listar[[i]] f0u[as1][[i]]  },{i,ni,lgr-2}];
t2=Table[{listar[[i]],-listar[[i]] (f0u[as2][[i]]  ) },{i,ni,lgr-2}];
t3=Table[{listar[[i]],-listar[[i]] (f0u[as3][[i]]  )},{i,ni,lgr-2}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" study f03"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]
 


ni=100;
t1=Table[{listar[[i]],-listar[[i]] (f1u[as1][[i]] ) },{i,ni,lgr-1}];
t2=Table[{listar[[i]],-listar[[i]] (f1u[as2][[i]]  ) },{i,ni,lgr-1}];
t3=Table[{listar[[i]],-listar[[i]] (f1u[as3][[i]]  )},{i,ni,lgr-1}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" study f13"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]
 


ni=100;
t1=Table[{listar[[i]],-listar[[i]] (f2u[as1][[i]] ) },{i,ni,lgr-1}];
t2=Table[{listar[[i]],-listar[[i]] (f2u[as2][[i]]  ) },{i,ni,lgr-1}];
t3=Table[{listar[[i]],-listar[[i]] (f2u[as3][[i]]   )},{i,ni,lgr-1}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" study f23"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]
 


ni=100;
t1=Table[{listar[[i]], listar[[i]]^2(wu[as1][[i]]) },{i,ni,lgr-1}];
t2=Table[{listar[[i]], listar[[i]]^2(wu[as2][[i]] ) },{i,ni,lgr-1}];
t3=Table[{listar[[i]], listar[[i]]^2(wu[as3][[i]] )},{i,ni,lgr-1}];


g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity,  Frame->True];

Print[" study w2"]
Show[g1,g2,g3,DisplayFunction->$DisplayFunction,Prolog->AbsolutePointSize[4]]
 







(* compute mass from asymptotics *)


ct=1;
cut=1;

ini=6; 

Do[

data=Table[{listar[[i]] ,f0u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{ 1/x ,1/x^2     } ,x];
cf01[k ]=Coefficient[u,1/x ];

data=Table[{listar[[i]] ,f1u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2       } ,x]; 
cf11[k ]=Coefficient[u,1/x ];


data=Table[{listar[[i]] ,f2u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2     } ,x];
cf21[k ]=Coefficient[u,1/x];


data=Table[{listar[[i]] ,wu[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x^2 ,1/x^3     } ,x];
cw2[k ]=Coefficient[u,1/x^2];

,{k,1,ny }]

f01=Table[cf01[k],{k,1,ny }] ;
f11=Table[cf11[k],{k,1,ny }] ;
f21=Table[cf21[k],{k,1,ny }] ;
w2=Table[cw2[k],{k,1,ny }] ;
 


1-Max[f01]/Min[f01]
1-Max[f11]/Min[f11]
1-Max[f21]/Min[f21]


f01
1+(2f01+f11 )/(f21)


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
 
(*
Series[-g[4,4],{r,Infinity,1}]
1+(2 const-rh)/r+O[1/r]^2
*)

const=Sum[f01[[i]],{i,1,ny}]/ny;

Mc=- const;

 
Mass= Mc;
Print["total Mass     = ",Mass];
 


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 crucial numerical test 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

test1=1+(2f01+f11 )/(f21);
ListPlot[test1,PlotRange->All, Frame->True,Axes->None]

err1=Max[Abs[test1]]


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA infinity -first derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


datx=ReadList["fx-inf.txt",{Number,Number ,Number ,Number,Number ,Number  }];
lung1=Length[datx] ;

r=Table[datx[[i]][[1]],{i,1,lung1}];
infF1x=Table[datx[[i]][[2]],{i,1,lung1}];
infF2x=Table[datx[[i]][[3]],{i,1,lung1}];
infF0x=Table[datx[[i]][[4]],{i,1,lung1}];
infZx=Table[datx[[i]][[5]],{i,1,lung1}];
infwx=Table[datx[[i]][[6]],{i,1,lung1}];


MassINF=Sum[infF0x[[i]],{i,1,ny}]/ny;
Print["Mass computed from the data at infinity=",MassINF];

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA infinity --second derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


datxx=ReadList["fxx-inf.txt",{Number,Number ,Number ,Number,Number ,Number  }];
lung1=Length[datxx] ;

r=Table[datxx[[i]][[1]],{i,1,lung1}];
 infWxx=Table[datxx[[i]][[6]],{i,1,lung1}];

JINF=Sum[ infWxx[[i]],{i,1,ny }]/(ny );
Print["J computed from the data at infinity=",JINF];
errJxx=1-Min[infWxx]/Max[infWxx];
Print["error Jinf=",errJxx];



infF1x/infF2x
infF1x/infF0x


Mass/MassINF//N


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA t-0 -- no conical singularities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


dat=ReadList["f-t0.txt",{Number,Number ,Number ,Number,Number ,Number  }];
lung1=Length[dat] ;

r=Table[dat[[i]][[1]],{i,1,lung1}];
t0F1=Table[dat[[i]][[2]],{i,1,lung1}];
t0F2=Table[dat[[i]][[3]],{i,1,lung1}];
t0F0=Table[dat[[i]][[4]],{i,1,lung1}];
t0W=Table[dat[[i]][[5]],{i,1,lung1}];
t0Z=Table[dat[[i]][[6]],{i,1,lung1}];

ratio= t0F1 -t0F2;
ListPlot[t0F1]
ListPlot[ratio]
Max[Abs[ratio]]



(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 computation MASS from the energy momentum tensor - Smarr relation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


asa1=2;
asa2=IntegerPart[ny/2];
asa3=ny-1;


 (* ordinea este: { Ttot,T44,T34}  *)
q=ReadList["T44.dat",{Number,Number,Number }];
lungq=Length[q];

diference=lungq-ntot;
Print["It must be zero! ",diference];
 

T34=Table[q[[k]][[3]],{k,1,lungq}]; 
ro=Table[q[[k]][[2]],{k,1,lungq}]; 
 Ttot=Table[q[[k]][[1]],{k,1,lungq}]; 

(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[ 
T34u[k]=Table[T34[[i]],{i,(k-1)*nx+1,k*nx}]; 
rou[k]=Table[ro[[i]],{i,(k-1)*nx+1,k*nx}]; 
Ttotu[k]=Table[Ttot[[i]],{i,(k-1)*nx+1,k*nx}]; 
(*	Print[T44u[k]];*)
,{k,1,ny}]

 
(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I plot the energy density for three different angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

ni=2;

t1=Table[{listar[[i]],-rou[asa1][[i]] },{i,ni,nx-7}];
t2=Table[{listar[[i]],-rou[asa2][[i]] },{i,ni,nx-7}];
t3=Table[{listar[[i]],-rou[asa3][[i]] },{i,ni,nx-7}];	
 
g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None ,PlotJoined->True];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ,PlotJoined->True, PlotStyle->{RGBColor[1,0,0]}];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity, Frame->True ,PlotJoined->True, PlotStyle->{RGBColor[1,1,0]}];
(* Print["black: theta=0 "];*)

 Print["profiles energy density"];
 Show[g1,g2,g3,DisplayFunction->$DisplayFunction];

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I plot the angular momentum density for three different angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
t1=Table[{listar[[i]],T34u[asa1][[i]] },{i,ni,nx-91}];
t2=Table[{listar[[i]],T34u[asa2][[i]] },{i,ni,nx-91}];
t3=Table[{listar[[i]],T34u[asa3][[i]] },{i,ni,nx-91}];
 
g1=ListPlot[t1,PlotRange->All,DisplayFunction->Identity, Frame->True,Axes->None ,PlotJoined->True,PlotRange->All];
g2=ListPlot[t2,PlotRange->All,DisplayFunction->Identity, Frame->True ,PlotJoined->True, PlotStyle->{RGBColor[1,0,0]}];
g3=ListPlot[t3,PlotRange->All,DisplayFunction->Identity, Frame->True ,PlotJoined->True, PlotStyle->{RGBColor[1,1,0]}];
(* Print["black: theta=0 "];*)

 Print["profiles T34"];
 Show[g1,g2,g3,DisplayFunction->$DisplayFunction];


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I compute the scalar field contribution to the total mass
& Smarr law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(* sqrt = E^(F0[r,t]+2 F1[r,t]+F2[r,t]) r Sqrt[g[r]] Sin[t]*)

(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)

Do[

 Mio2[k]=Table[{listar[[i]], E^(f0u[k][[i]]+2 f1u[k][[i]]+f2u[k][[i]])  listar[[i]]^2  Ttotu[k][[i]]},{i,ni,nx-1}];

Jio2[k]=Table[{listar[[i]], E^(f0u[k][[i]]+2 f1u[k][[i]]+f2u[k][[i]])  listar[[i]]^2  T34u[k][[i]]},{i,ni,nx-1}];


,{k,2,ny-1}];



(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)
Do[
Ma2[k]=ListIntegrate[Mio2[k],2]//N;
Ja2[k]=ListIntegrate[Jio2[k],2]//N;
 ,{k,2,ny-1}];

 
Ma2[1]=Ma2[2];
Ma2[ny]=Ma2[ny-1];
 
Ja2[1]=Ja2[2];
Ja2[ny]=Ja2[ny-1];


 Minn2=Table[{unghi[[k]], Sin[unghi[[k]]] Ma2[k]},{k,1,ny}];
Jinn2=Table[{unghi[[k]], Sin[unghi[[k]]] Ja2[k]},{k,1,ny}];

 
 Mint=ListIntegrate[Minn2,2]//N; 
Jint=ListIntegrate[Jinn2,2]//N; 
 Print["Mintegral = ",Mint];  
 Print["Jintegral = ",Jint]; 

(*max, r(max)*)
d1=ReadList["sup.dat",{Number,Number    }]
Zm=d1[[1]][[1]];
rm= d1[[1]][[2]];




- MassINF/Mint/alfa^2


 JINF/Jint/4/alfa^2






asa= Table[{ alfa,wf,c1,c2,c3,MassINF,JINF,Mint,Jint,minf0,f0H,f1H,maxW ,Zm,rm}] 








stmp=OpenAppend["tmp.txt"];
Write[stmp,asa];
Close[stmp] ;



