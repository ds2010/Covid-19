$TITLE  Covid-19 application: First wave

* This script is developed for non-commercial use only. All rights reserved.

* Timo Kuosmanen, Aaron Tan, and Sheng Dai,
* Performance of English NHS hospitals during the first and second waves of the COVID-19 pandemic (February 5, 2021).
* Available at Researchgate:


SETS
     i       "DMU's"   /i1*i1938/
     j       'inputs and outputs' /y, x1, x2, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11/
     inp(j)  'inputs'  /x1, x2/
     outp(j) 'outputs' /y/
     zvar(j) 'z variables' /z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11/
;

alias(i,h)
;

# Import all data
Parameter data(i,j) ;
$GDXIN covid_1_18.gdx
$LOAD data
$GDXIN
;

PARAMETERS
X(i,inp)       'inputs of firm i'
Y(i)           'outputs of firm i'
Z(i, zvar)     "Environmental indicator"
Evalue(i)      "residual"
Bvalue(i,inp)  'output beta values'
Dvalue(zvar)   "Estimated Delta-coefficient"
Avalue(i)
;
X(i,inp) = data(i,inp);
Y(i)     = data(i,'y');
Z(i,zvar) = data(i,zvar);

VARIABLES
D(zvar)         "Coefficient of the z-variable"
E(i)            "resduals"
SSE             "Sum of squares of residuals"
alpha(i)
;

POSITIVE VARIABLES
B(i,inp)        beta coefficients
Yhat(i)
;

**** formula's for calculating CNLS***
EQUATIONS
QSSE             objective=sum of squares of errors
QREGRESSION(i)   regression equation
Qlog(i)
QCONC(i,h)
beta1(i)
beta2(i)
;

QSSE..           SSE=e=sum(i,E(i)*E(i));
QREGRESSION(i).. log(Y(i)) =e= log(Yhat(i)+1) + sum(zvar, D(zvar)*Z(i,zvar)) + E(i);
Qlog(i)..        Yhat(i) =e= alpha(i) + sum(inp, B(i,inp)*X(i,inp)) - 1;
QCONC(i,h)..     alpha(i) + sum(inp, B(i,inp)*X(i,inp)) =l= alpha(h) + sum(inp, B(h,inp)*X(i,inp));
beta1(i)..       B(i,'x1') =l= 1;
beta2(i)..       B(i,'x2') =l= 1;

MODEL covid /all/;

*solver configuration
OPTION solvelink = 0;
OPTION limrow    = 0;
OPTION limcol    = 0;
OPTION SOLPRINT  = OFF;
OPTION optcr     = 0.0;
OPTION iterlim   = 10000000;
OPTION reslim    = 10000000;
OPTION decimals=8;
Option NLP = KNITRO;

SOLVE covid using NLP Minimizing SSE;

Evalue(i)=E.l(i);
Bvalue(i,inp)=B.l(i,inp);
Dvalue(zvar)=D.l(zvar);
Avalue(i)= alpha.l(i);

Execute_unload "vrs_covid_1.gdx" Evalue Bvalue Dvalue Avalue
