$TITLE  Covid-19 application: convex expectile regression (CER) by CNLS+ Algorithm

* This script is developed for non-commercial use only. All rights reserved.

* Timo Kuosmanen, Aaron Tan, and Sheng Dai,
* Performance of English NHS hospitals during the first and second waves of the COVID-19 pandemic (February 5, 2021).
* Available at Researchgate:


SETS
i                observations /i1*i3260/
k                outputs /k1*k2/
m                index /m1*m1/
l                zvars /l1*l11/
;

alias(i,j)
;

Scalars
Activetmp        control the loop in CNLS+G
Activetmp1       control the sub-loop in CNLS+G
tau              expectile   /0.05/
;

* Output
PARAMETERS   C(i,m);
$GDXIN Yvar.gdx
$LOAD C
$GDXIN
;
* Inputs
PARAMETERS   Y(i,k);
$GDXIN Xvar.gdx
$LOAD Y
$GDXIN
;
* Z variables
PARAMETERS   Z(i,l) ;
$GDXIN Zvar.gdx
$LOAD Z
$GDXIN
;

PARAMETERS
Active(i,j)      active (added) violated concavity constraint by iterative procedure
Active2(i,j)     violated concavity constraint
Cutactive(i,j)   active (added) concavity constraint by sweet spot
distance(i,j)    distance matrix for sweet spot
distcut(i,m)     3rd percentile of the distances by sweet spot
;

$GDXIN dist.gdx
$LOAD distance
$GDXIN
;
$GDXIN discut.gdx
$LOAD distcut
$GDXIN
;

VARIABLES
SS1              sum of square of errors
D1(l)
SS2              sum of square of errors
D2(l)
A1(i)
A2(i)
;

POSITIVE VARIABLES
B1(i,k)           beta coefficients (positivity quarantees monotonicity)
Chat1(i)
B2(i,k)           beta coefficients (positivity quarantees monotonicity)
Chat2(i)

eplus1(i)         error term  
eminus1(i)
eplus2(i)
eminus2(i)
;

**** formula's for calculating CSLS***
EQUATIONS
QSSE1             objective=sum of squares of errors
QREGRESSION1(i)   regression equation
Qlog1(i)
QCONC1(i)
QCONB1(i,j)
beta11(i)
beta12(i)
;

QSSE1..                          SS1=e=(1-tau)*sum(i,eminus1(i)*eminus1(i)) + tau*sum(i,eplus1(i)*eplus1(i));
QREGRESSION1(i)..                log(C(i,'m1')) =e= log(Chat1(i) + 1) + sum(l, D1(l)*Z(i,l)) - eminus1(i) + eplus1(i);
Qlog1(i)..                       Chat1(i) =e= sum(k, B1(i,k)*Y(i,k)) + A1(i) - 1;
QCONC1(i)..                      sum(k, B1(i,k)*Y(i,k))+ A1(i) =l= sum(k, B1(i++1,k)*Y(i,k)) + A1(i++1);
QCONB1(i,j)$(Cutactive(i,j))..   sum(k, B1(i,k)*Y(i,k))+ A1(i) =l= sum(k, B1(j,k)*Y(i,k)) + A1(j);
beta11(i)..       		         B1(i,'k1') =l= 1;
beta12(i)..       		         B1(i,'k2') =l= 1;
Model CNLS1 /QSSE1, QREGRESSION1, Qlog1, QCONC1, QCONB1, beta11, beta12/

EQUATIONS
QSSE2             objective=sum of squares of errors
QREGRESSION2(i)   regression equation
Qlog2(i)
QCONC2(i)
QCONC3(i,j)
QCONB2(i,j)
beta21(i)
beta22(i)
;

QSSE2..                          SS2 =e= (1-tau)*sum(i,eminus2(i)*eminus2(i)) + tau*sum(i,eplus2(i)*eplus2(i));
QREGRESSION2(i)..                log(C(i,'m1')) =e= log(Chat2(i) + 1) + sum(l, D2(l)*Z(i,l)) - eminus2(i) + eplus2(i);
Qlog2(i)..                       Chat2(i) =e= sum(k, B2(i,k)*Y(i,k)) + A2(i) - 1;
QCONC2(i)..                      sum(k, B2(i,k)*Y(i,k)) + A2(i) =l= sum(k, B2(i++1,k)*Y(i,k)) + A2(i++1);
QCONC3(i,j)$(Cutactive(i,j))..   sum(k, B2(i,k)*Y(i,k)) + A2(i) =l= sum(k, B2(j,k)*Y(i,k)) + A2(j);
QCONB2(i,j)$(Active(i,j))..      sum(k, B2(i,k)*Y(i,k)) + A2(i) =l= sum(k, B2(j,k)*Y(i,k)) + A2(j);
beta21(i)..       		         B2(i,'k1') =l= 1;
beta22(i)..       		 	     B2(i,'k2') =l= 1;
Model CNLS2 /QSSE2, QREGRESSION2, Qlog2, QCONC2, QCONC3, QCONB2, beta21, beta22/
;

*solver configuration
CNLS2.workfactor = 1.5;
OPTION solvelink = 0;
OPTION limrow    = 0;
OPTION limcol    = 0;
OPTION SOLPRINT  = OFF;
OPTION optcr     = 0.0;
OPTION iterlim   = 10000000;
OPTION reslim    = 10000000;
OPTION decimals=8;
OPTION NLP=KNITRO;

Cutactive(i,j)=yes$(distance(i,j)<= distcut(i,'m1'));

SOLVE CNLS1 using NLP Minimizing SS1;

loop((i,j),
     Active(i,j)=0;
);

*go into the loop
Activetmp1=0;
loop(i,
     Activetmp=0;
     loop(j,
          Active2(i,j) = sum(k, B1.l(i,k)*Y(i,k)) + A1.l(i) - sum(k, B1.l(j,k)*Y(i,k)) - A1.l(j);
          Activetmp$(Active2(i,j) GT Activetmp)=Active2(i,j);
     );
*find the maximal violated constraint in sub-loop and added into the active matrix
     Active(i,j)=yes$(Active2(i,j)>= Activetmp and Activetmp >0);
     Activetmp1$(Activetmp GT Activetmp1)=Activetmp;
);

while(Activetmp>0.0001,

SOLVE CNLS2 using NLP Minimizing SS2;

Activetmp1=0;
*go into the loop
loop(i,
     Activetmp=0;
*go into the sub-loop and find the violated concavity constraints
     loop(j,
          Active2(i,j) = sum(k, B2.l(i,k)*Y(i,k)) + A2.l(i) - sum(k, B2.l(j,k)*Y(i,k)) - A2.l(j);
          Activetmp$(Active2(i,j) GT Activetmp)=Active2(i,j);
     );
*find the maximal violated constraint in sub-loop and added into the active matrix
     Active(i,j)=yes$(Active(i,j))+yes$(active2(i,j)>= Activetmp and Activetmp >0);
     Activetmp1$(Activetmp GT Activetmp1)=Activetmp;
     );
);

* Export estimation results
Execute_unload "Results05.gdx" A2.l, B2.l, D2.l, eplus2.l, eminus2.l
