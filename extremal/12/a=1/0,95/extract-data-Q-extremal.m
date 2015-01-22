(* ::Package:: *)

(* 8 Martie 14 *)


 (*<<NumericalMath`ListIntegrate` *)


Off[SetDelayed::"write"]
<<NumericalMath`ListIntegrate`


(* *)


Off[General::spell1]
Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];

(* 1  2   3   4  5  6   7*)
(* nr,w,alpha,c1,c2,c3,rh*)
conf=ReadList["res.txt",{Number,Number ,Number ,Number ,Number,Number,Number }]

 
nr=conf[[1]][[1]];
rh= conf[[1]][[7]];
alfa= conf[[1]][[3]];
c1= conf[[1]][[4]];
c2= conf[[1]][[5]];
c3= conf[[1]][[6]]; 
w= conf[[1]][[2]]
 

Print["winding number n = ",nr];
Print["rh   = ",rh];
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

a=ReadList["functf.dat",{Number,Number,Number,Number ,Number   }];


lung1=Length[a];


(*Datele sunt salvate direct cu indexarea globala*)
F1=Table[a[[k]][[1]],{k,1,lung1}];
F2=Table[a[[k]][[2]],{k,1,lung1}];
F0=Table[a[[k]][[3]],{k,1,lung1}];  
Z=Table[a[[k]][[4]],{k,1,lung1}];
W=Table[a[[k]][[5]],{k,1,lung1}]; 
  
(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[

F1u[k]=Table[F1[[i]],{i,(k-1)*nx+1,k*nx}];
F2u[k]=Table[F2[[i]],{i,(k-1)*nx+1,k*nx}];
F0u[k]=Table[F0[[i]],{i,(k-1)*nx+1,k*nx}];  
Zu[k]=Table[Z[[i]],{i,(k-1)*nx+1,k*nx}]; 
Wu[k]=Table[W[[i]],{i,(k-1)*nx+1,k*nx}]; 

,{k,1,ny}];
 
 
as1=2;
as2=IntegerPart[ny/2];
as3=ny-1;

sa1=3;
sa2=IntegerPart[nx/2];
sa3=nx-1;


Print["rmax = ",gr[[nx]][[1]]];


listar


minF0=Min[F0] 
maxF0=Max[ F0 ] 
minF1=Min[F1] 
maxF1=Max[ F1 ] 
minF2=Min[F2] 
maxF2=Max[ F2 ] 
minW=Min[W] 
maxW=Max[ W ] 
maxZ=Max[ Z ] 
minZ=Min[Z] 

f0H=F0u[1][[1]]
f1H=F1u[1][[1]]
f2H=F2u[1][[1]]
ZH=Zu[ny][[1]]

 



n1=nx-123;

t1=Table[ {listar[[i]] ,E^(2*F0u[as1][[i]])},{i,2,n1}];
t2=Table[{ listar[[i]]  ,E^(2*F0u[as2][[i]])},{i,2,n1}];
t3=Table[{ listar[[i]] ,E^(2*F0u[as3][[i]])},{i,2,n1}];



as=Table[{  F0u[as3][[i]]},{i,2,n1}] ;
Min[F0]


E^( 2Min[F0]) 


cut=2;

t1=Table[{i,F1u[as1][[i]]},{i,1,lgr}];
t2=Table[{i,F1u[as2][[i]]},{i,1,lgr}];
t3=Table[{i,F1u[as3][[i]]},{i,1,lgr}];



Print[" FUNCTION F1"]

t1=Table[ {listar[[i]] ,F1u[as1][[i]]},{i,1,lgr-cut}];
t2=Table[{ listar[[i]]  ,F1u[as2][[i]]},{i,1,lgr-cut}];
t3=Table[{ listar[[i]] ,F1u[as3][[i]]},{i,1,lgr-cut}];





   


t1=Table[{i,F2u[as1][[i]]},{i,1,lgr}];
t2=Table[{i,F2u[as2][[i]]},{i,1,lgr}];
t3=Table[{i,F2u[as3][[i]]},{i,1,lgr}];



Print[" FUNCTION F2"]


t1=Table[ {listar[[i]] ,F2u[as1][[i]]},{i,1,lgr-cut}];
t2=Table[{ listar[[i]]  ,F2u[as2][[i]]},{i,1,lgr-cut}];
t3=Table[{ listar[[i]] ,F2u[as3][[i]]},{i,1,lgr-cut}];




t1=Table[{i,F0u[as1][[i]]},{i,1,lgr}];
t2=Table[{i,F0u[as2][[i]]},{i,1,lgr}];
t3=Table[{i,F0u[as3][[i]]},{i,1,lgr}];



Print[" FUNCTION F0"]


t1=Table[ {listar[[i]] ,F0u[as1][[i]]},{i,2,lgr-cut}];
t2=Table[{ listar[[i]]  ,F0u[as2][[i]]},{i,2,lgr-cut}];
t3=Table[{ listar[[i]] ,F0u[as3][[i]]},{i,2,lgr-cut}];




t1=Table[{i,Zu[as1][[i]] },{i,1,lgr}];
t2=Table[{i,Zu[as2][[i]] },{i,1,lgr}];
t3=Table[{i,Zu[as3][[i]] },{i,1,lgr}];



Print[" FUNCTION Z"]

t1=Table[ {listar[[i]] ,Zu[as1][[i]]},{i,2,lgr-cut}];
t2=Table[{ listar[[i]]  ,Zu[as2][[i]]},{i,2,lgr-cut}];
t3=Table[{ listar[[i]] ,Zu[as3][[i]]},{i,2,lgr-cut}];


 


t1=Table[{i,Wu[as1][[i]] },{i,1,lgr}];
t2=Table[{i,Wu[as2][[i]] },{i,1,lgr}];
t3=Table[{i,Wu[as3][[i]] },{i,1,lgr}];



Print[" FUNCTION W"]

t1=Table[ {listar[[i]] ,Wu[as1][[i]]},{i,1,lgr-cut}];
t2=Table[{ listar[[i]]  ,Wu[as2][[i]]},{i,1,lgr-cut}];
t3=Table[{ listar[[i]] ,Wu[as3][[i]]},{i,1,lgr-cut}];;


 







cut1= 2;

t1=Table[ {listar[[i]] ,Wu[as1][[i]] listar[[i]]^1},{i,1,lgr-cut1}];
t2=Table[{ listar[[i]]  ,Wu[as2][[i]] listar[[i]]^1},{i,1,lgr-cut1}];
t3=Table[{ listar[[i]] ,Wu[as3][[i]] listar[[i]]^1},{i,1,lgr-cut1}];;




nr1=170;
t1=Table[ {listar[[i]] ,-Zu[as3][[i]]},{i,2,nr1}];
t2=Table[{ listar[[i]]  ,F0u[as3][[i]]},{i,2,nr1}]; 





(*Do[
Print["z =",i,"  ",unghi[[i]]];
,{i,1,ny}]*)


ni=200;

cut= 1;
t1=Table[{listar[[i]], listar[[i]] (F0u[as1][[i]] ) },{i,ni,lgr- cut}];
t2=Table[{listar[[i]], listar[[i]] (F0u[as2][[i]]  ) },{i,ni,lgr- cut}];
t3=Table[{listar[[i]], listar[[i]] (F0u[as3][[i]]   )},{i,ni,lgr- cut}];



Print[" FUNCTION F01(theta)"]
 


 t1=Table[{listar[[i]], listar[[i]] (F1u[as1][[i]] ) },{i,ni,lgr-cut}];
t2=Table[{listar[[i]], listar[[i]] (F1u[as2][[i]]  ) },{i,ni,lgr-cut}];
t3=Table[{listar[[i]], listar[[i]] (F1u[as3][[i]]   )},{i,ni,lgr-cut}];



Print[" FUNCTION F11(theta)"]
 


cut


 t1=Table[{listar[[i]], listar[[i]] (F2u[as1][[i]] ) },{i,ni,lgr-cut}];
t2=Table[{listar[[i]], listar[[i]] (F2u[as2][[i]]  ) },{i,ni,lgr-cut}];
t3=Table[{listar[[i]], listar[[i]] (F2u[as3][[i]]   )},{i,ni,lgr-cut}];



Print[" FUNCTION F21(theta)"]
 


 t1=Table[{listar[[i]], listar[[i]]^1 (Wu[as1][[i]] ) },{i,ni,lgr-cut}];
t2=Table[{listar[[i]], listar[[i]]^1 (Wu[as2][[i]]  ) },{i,ni,lgr-cut}];
t3=Table[{listar[[i]], listar[[i]]^1 (Wu[as3][[i]]   )},{i,ni,lgr-cut}];



Print[" FUNCTION F21(theta)"]
 


nf=nx- cut;
ti=Table[{i, listar[[nf]] (F0u[i][[nf]] ) },{i,1,ny}]


nf=nx- cut;
ti=Table[{i, listar[[nf]]^1 (Wu[i][[nf]] ) },{i,1,ny}]


ct=1;
cut=1;

ini=5; 

Do[

data=Table[{listar[[i]] ,F0u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{ 1/x ,1/x^2     } ,x];
cf01[k ]=Coefficient[u,1/x ];

data=Table[{listar[[i]] ,F1u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2       } ,x]; 
cf11[k ]=Coefficient[u,1/x ];


data=Table[{listar[[i]] ,F2u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2     } ,x];
cf21[k ]=Coefficient[u,1/x];

 data=Table[{listar[[i]] ,Wu[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{ 1/x ,1/x^2     } ,x];
cW[k ]=Coefficient[u,1/x];

,{k,1,ny }]

f01=Table[cf01[k],{k,1,ny }] ;
f11=Table[cf11[k],{k,1,ny }] ;
f21=Table[cf21[k],{k,1,ny }] ;
W3=Table[cW[k],{k,1,ny }] ;

 constJINF=Sum[W3[[i]],{i,1,ny}]/ny


1-Max[f01]/Min[f01]
1-Max[f11]/Min[f11]
1-Max[f21]/Min[f21]
1-Max[W3]/Min[W3]


f01
1+(2f01+f11 )/(f21)






(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 crucial numerical test 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

test1=1+(2f01+f11 )/(f21);

err1=Max[Abs[test1]]


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA infinity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


dat=ReadList["fx-inf.txt",{Number,Number ,Number ,Number,Number ,Number  }];
lung1=Length[dat] ;

r=Table[dat[[i]][[1]],{i,1,lung1}];
infF1=Table[dat[[i]][[2]],{i,1,lung1}]
infF2=Table[dat[[i]][[3]],{i,1,lung1}];
infF0=Table[dat[[i]][[4]],{i,1,lung1}]
infZ=Table[dat[[i]][[5]],{i,1,lung1}];



w


constINF=Sum[infF0[[i]],{i,1,ny}]/ny










(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
computation Mass from asymptotics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(* th result:*) 
(*
Series[-g[4,4],{r,Infinity,1}]
1+(2 const-rh)/r+O[1/r]^2: non-extremal

Series[gtt,{r,Infinity,1}]
-1+(-2 ct+2 rh)/r+O[1/r]^2: extremal
*)
(* atentie!! la extremal este diferit! *)


const=Sum[f01[[i]],{i,1,ny}]/ny;

Mc=  constINF;

(* non-extremal:
Mc=  constINF;
MSch=rh/2 ;
Mass=MSch+Mc;
*)


MSch=rh ;
Mass=MSch+Mc;

Print["Mass Schw     = ",MSch];
Print["Mass correction= ",Mc];
Print["total Mass     = ",Mass];

(*Print[ MSAdS/Mc ]*)


(* difference mass computed at infinity/mass interpolated *)
1-Abs[const/constINF] 






(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA t-0 -- no conical singularities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


dat=ReadList["f-t0.txt",{Number,Number ,Number ,Number,Number,Number  }];
lung1=Length[dat] 

r=Table[dat[[i]][[1]],{i,1,lung1}];
t0F1=Table[dat[[i]][[2]],{i,1,lung1}];
t0F2=Table[dat[[i]][[3]],{i,1,lung1}];
t0F0=Table[dat[[i]][[4]],{i,1,lung1}];
t0Z=Table[dat[[i]][[5]],{i,1,lung1}];

ratio= t0F2-t0F1;
Max[ratio]




(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE HORIZON DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


dat=ReadList["f-0.txt",{Number,Number ,Number ,Number,Number ,Number  }];
lung1=Length[dat] ;

unghi=Table[dat[[i]][[1]],{i,1,lung1}];
hF1=Table[dat[[i]][[2]],{i,1,lung1}]; 
hF2=Table[dat[[i]][[3]],{i,1,lung1}]; 
hF0=Table[dat[[i]][[4]],{i,1,lung1}]; 
hZ=Table[dat[[i]][[5]],{i,1,lung1}]
 




(* angular dependence of the entropy corrections *)






(*%%%%%%% event horizon area %%%%%%%%*)
AH0=4 Pi rh^2;

iAHc=Table[{unghi[[k]],1/2  Sin[unghi[[k]]] E^((hF1[[k]] +hF2[[k]] )) },{k,1,lung1 }];

(* 2 because I integrate between 0, Pi/2 *)
 AHc= 2Integrate[Interpolation[iAHc,InterpolationOrder->1][x],{x,0,Pi/2}];

AH=AH0 AHc;

Print["Schw area  AH0= ",AH0];
Print["correction AHc= ",AHc];
Print["Event horizon area =",AH];







(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
computation Le, Lp -- see MATH code for derivation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(*Le=2 E^F2[0,\[Pi]/2] \[Pi] rH*)
Le=2 E^hF2[[ny]] \[Pi] rh;
Print["Le = ",Le ];

(*Lp=2 \!\(
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Pi]\)]\(\(
\*SuperscriptBox[\(E\), \(F1[0, t]\)]\ rH\) \[DifferentialD]t\)\)*)
Lp1= Table[{unghi[[k]],2 rh    E^hF1[[k]]   },{k,1,lung1 }];
Lp=  2ListIntegrate[Lp1,2]//N;(*factor 2: because I integrate [0,pi/2]*) 
Print["Lp = ",Lp ];
Print["  " ];








(* computation MASS from the energy momentum tensor - Smarr relation *)

asa1=2;
asa2=IntegerPart[ny/2];
asa3=ny-1;


 (* ordinea este: { T34,T44,Ttot }  *)
q=ReadList["T44.dat",{Number,Number,Number }];
lungq=Length[q];

diference=lungq-ntot;
Print["It must be zero! ",diference];
 

T34=Table[q[[k]][[1]],{k,1,lungq}]; 
ro=Table[q[[k]][[2]],{k,1,lungq}]; 
 Ttot=Table[q[[k]][[3]],{k,1,lungq}]; 

(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[ 
 
(* T34u[k]=Table[T34[[i]],{i,(k-1)*nx+1,k*nx}]; *)
rou[k]=Table[ro[[i]],{i,(k-1)*nx+1,k*nx}]; 
Ttotu[k]=Table[Ttot[[i]],{i,(k-1)*nx+1,k*nx}]; 
(*	Print[T44u[k]];*)
,{k,1,ny}]


(* I redefine the angular momentum density *)
eps1=0.0001;

Do[ 
(* H[r_]:=rh+Sqrt[rh^2+r^2];*)
(* generic expression theory: T34u[k]=-((2 E^(-2 F0[r,t]) nr H[r]^2 (w g[r]+nr W1[r,t]) Z[r,t]^2)/r^4);*)
  T34u[k]=Table[-(2 E^(-2F0u[k][[i]]) nr (rh+Sqrt[rh^2+listar[[i]]^2])^2 ( w (listar[[i]]^2+rh^2)+nr Wu[k][[i]]) Zu[k][[i]]^2)/(listar[[i]]+eps1)^4,{i,1,nx}];
 
 (*	Print[T44u[k]];*)
,{k,1,ny}]



(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I plot the energy density for three different angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

ni=2;
cut=143;

t1=Table[{listar[[i]],-rou[asa1][[i]] },{i,ni,nx-cut}];
t2=Table[{listar[[i]],-rou[asa2][[i]] },{i,ni,nx-cut}];
t3=Table[{listar[[i]],-rou[asa3][[i]] },{i,ni,nx-cut}];
 
(* Print["black: theta=0 "];*)

 Print["profiles energy density"];


t1=Table[{i,-rou[as1][[i]] },{i,1,lgr}];
t2=Table[{i,-rou[as2][[i]] },{i,1,lgr}];
t3=Table[{i,-rou[as3][[i]] },{i,1,lgr}];



Print[" energy density"]

t1=Table[ {listar[[i]] ,-rou[as1][[i]]},{i,1,lgr-cut}];
t2=Table[{ listar[[i]]  ,-rou[as2][[i]]},{i,1,lgr-cut}];
t3=Table[{ listar[[i]] ,-rou[as3][[i]]},{i,1,lgr-cut}];;


 


ni



cut=115;


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I plot the constraint T34 for three different angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
t1=Table[{listar[[i]],T34u[asa1][[i]] },{i,ni,nx-cut}];
t2=Table[{listar[[i]],T34u[asa2][[i]] },{i,ni,nx-cut}];
t3=Table[{listar[[i]],T34u[asa3][[i]] },{i,ni,nx-cut}];	
 
(* Print["black: theta=0 "];*)

 Print["profiles T34"];





(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I compute the scalar field contribution to the total mass
& Smarr law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(* sqrt = E^(F0[r,t]+2 F1[r,t]+F2[r,t]) r Sqrt[g[r]] Sin[t]*)

(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)

Do[

 Mio2[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2]  Ttotu[k][[i]]},{i,ni,nx-1}];

 Mio3[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2]  T34u[k][[i]]},{i,ni,nx-1}];

,{k,2,ny-1}];



(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)
Do[
Ma2[k]=ListIntegrate[Mio2[k],2]//N;
Ma3[k]=ListIntegrate[Mio3[k],2]//N;
 ,{k,2,ny-1}];

 
Ma2[1]=Ma2[2];
Ma2[ny]=Ma2[ny-1];
 
Ma3[1]=Ma3[2];
Ma3[ny]=Ma3[ny-1];

 Minn2=Table[{unghi[[k]], Sin[unghi[[k]]] Ma2[k]},{k,1,ny}];
  Minn3=Table[{unghi[[k]], Sin[unghi[[k]]] Ma3[k]},{k,1,ny}];

 Mint=ListIntegrate[Minn2,2]//N; 
Jint=ListIntegrate[Minn3,2]//N; 
 Print["Mintegral = ",Mint];  
 Print["Jintegral = ",Jint];  



 constJINF/2/Jint


(* Smarr law Schw BH;
MSch-2 TH0 1/4 AH0 
 *)

Print["total mass     = ",Mass];
Print["total J     = ",1/2 constJINF ];
Print["mass integral  = ",-Mint 1/2 2alfa^2];
Print["J integral  = ",Jint];
Print["AH    = ",  AH];
Print["factor TH S    = ",2 TH 1/4 AH];
Print["factor OmegaHJ = ",2  w   (1/2 constJINF -Jint)];

Print[" "];

(* smarr relation *)
errSmarr=1-(Mass +2  w   (1/2 constJINF -Jint) )/(-Mint 1/2 2alfa^2);
Print[" errSmarr= ",errSmarr]; 






(* difference mass computed at infinity/mass interpolated *)
1-Abs[const/constINF] 





asa= Table[{rh,alfa,c1,c2,c3,Jint,TH,Mass,AH,err1,minF0,f0H,f1H,f2H,ZH,Mint,Le,Lp,w,constJINF,maxF0,minF1,maxF1,minF2,maxF2,minW,maxW,minZ,maxZ,errSmarr }] 





stmp=OpenAppend["tmp.txt"];
Write[stmp,asa];
Close[stmp] ;
 
