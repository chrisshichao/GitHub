/*************   Cointegration Analysis of Real Money Balance     *************/
/*************                      Chao Shi                      *************/

/******************************************************************************/

clear				/* drop all data and value labels from memory */
clear matrix
set memory 100m		/* reset the amount of memory allocated to Stata */ 


			/*** 0. Import Data ***/

/* A. Import data in .csv format */
insheet using "C:\Users\chris\Desktop\ECON825\Money.csv", names comma clear

/* missing values in .csv should be coded as ".", not "n.a." or " " */
/* record all variable names in lower case b/c Stata is case sensitive */
describe		/* check numeric vs. string variables */

/* Declare the data to be time series & specify the time variable */    
cap drop time
gen time=quarterly(date,"yq")	
tsset time, quarterly	
describe		
list date time mp y r	/* check missing values: none */

/* B. Label log transformed variables */
label var mp "log (real money balance M2)"
label var y "log (real private output)"
label var r "interest rate"
codebook mp y r

/******************************************************************************/

			/*** 1. Unit Root Tests ***/

	/* 1.1 MP: log of real money balance M2 */

/* A. time series plot for level and first difference */
tsline mp
tsline d.mp

gen trend=_n   /* this step is necessary b/c _n cannot be used directly as a regressor */

/* B1. ADF Test for mp */	

/* Schwert (1989), pmax=int[12*{(T+1)/100}^0.25]=13 */
scalar T=168
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=0
scalar diff=0
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

reg mp trend l.mp l(1/13)d.mp 
estat ic
matrix mp_order=r(S)
reg mp trend l.mp l(1/12)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/11)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/10)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/9)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/8)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/7)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/6)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/5)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/4)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/3)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/2)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp l(1/1)d.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg mp trend l.mp  if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))

matlist mp_order   /* AIC: P=1 */
dfuller mp, trend lag(1) /* z=-3.161  5% CV=-3.441 ==> failed to reject H0: mp is DS with drift */

/* B2. ADF Test for d.mp */	

/* Schwert (1989), pmax=int[12*{(T+1)/100}^0.25]=13 */
scalar T=168
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=0
scalar diff=1
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

reg d.mp ld.mp l(1/13)d2.mp 
estat ic
matrix mp_order=r(S)
reg d.mp ld.mp l(1/12)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/11)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/10)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/9)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/8)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/7)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/6)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/5)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/4)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/3)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/2)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp l(1/1)d2.mp if _n>p_reg
estat ic
matrix mp_order=(mp_order \ r(S))
reg d.mp ld.mp  if _n>p_reg 
estat ic
matrix mp_order=(mp_order \ r(S))

matlist mp_order   /* AIC: P=0 */
dfuller d.mp

/* C1. ADF-GLS Tests for mp */
dfgls mp				/* AIC: P=1 */		
dfgls mp, maxlag(1)

/* C2. ADF-GLS Tests for d.mp */
dfgls d.mp, notrend	/* AIC: P=4 */
dfgls d.mp, maxlag(4) notrend


	/* 1.2 Y: log of real private output */

/* A. time series plot for level and first difference */
tsline y
tsline d.y

/* B1. ADF Test for y */	

scalar T=168
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=0
scalar diff=0
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

reg y trend l.y l(1/13)d.y if _n>p_reg
estat ic
matrix y_order=r(S)
reg y trend l.y l(1/12)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/11)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/10)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/9)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/8)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/7)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/6)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/5)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/4)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/3)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/2)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y l(1/1)d.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg y trend l.y if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))

matlist y_order /* AIC: P=2 */
dfuller y, trend lag(2)

/* C1. ADF-GLS Tests for y */
dfgls y	
dfgls y, maxlag(13)

/* C2. ADF-GLS Tests for d.y */
dfgls d.y, notrend
dfgls d.y, maxlag(1) notrend


	/* 1.3 R: interest rate */

/* A. time series plot for level and first difference */
tsline r
tsline d.r

/* B1. ADF Test for r */
scalar T=168
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=0
scalar diff=0
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg


reg r l.r l(1/13)d.r if _n>p_reg
estat ic
matrix r_order=r(S)
reg r l.r l(1/12)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/11)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/10)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/9)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/8)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/7)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/6)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/5)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/4)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/3)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/2)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r l(1/1)d.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))
reg r l.r if _n>p_reg
estat ic
matrix r_order=(r_order \ r(S))		

matlist r_order /* AIC: P=7*/
dfuller r, lag(7) 

/* B2. ADF Test for d.r */
scalar T=168
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=0
scalar diff=1
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

reg d.r ld.r l(1/13)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=r(S)
reg d.r ld.r l(1/12)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/11)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/10)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/9)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/8)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/7)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/6)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/5)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/4)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/3)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/2)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r l(1/1)d2.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))
reg d.r ld.r if _n>p_reg, noconstant
estat ic
matrix r_order=(r_order \ r(S))

matlist r_order /* AIC: P=6 */
dfuller d.r, lag(6) noconstant

/* C1. ADF-GLS Tests for r */
dfgls r	/* AIC: P=7*/
dfgls r, maxlag(7)

/******************************************************************************/

			/*** 2. Residual-Based Cointegration Test ***/

scalar T=168
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=0
scalar diff=0		
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

	/* 2.1 MP as dependent variable */
	
reg mp y r
cap drop res
predict res, residuals

reg res l.res l(1/13)d.res, noconst
estat ic
matrix res_order=r(S)
reg res l.res l(1/12)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/11)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/10)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/9)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/8)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/7)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/6)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/5)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/4)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/3)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/2)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res ld.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))

matlist res_order				/* AIC: P=1*/
dfuller res, lag(1) noconstant

	/* 2.2 MP as dependent variable */
reg y mp r
cap drop res
predict res, residuals

reg res l.res l(1/13)d.res, noconst
estat ic
matrix res_order=r(S)
reg res l.res l(1/12)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/11)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/10)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/9)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/8)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/7)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/6)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/5)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/4)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/3)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/2)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res ld.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))

matlist res_order				/* AIC: P=9*/
dfuller res, lag(9) noconstant

	/* 2.2 R as dependent variable */
reg r mp y
cap drop res
predict res, residuals

reg res l.res l(1/13)d.res, noconst
estat ic
matrix res_order=r(S)
reg res l.res l(1/12)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/11)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/10)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/9)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/8)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/7)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/6)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/5)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/4)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/3)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res l(1/2)d.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res ld.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))
reg res l.res if _n>p_reg, noconst
estat ic
matrix res_order=(res_order \ r(S))

matlist res_order				/* AIC: P=8 */
dfuller res, lag(8) noconstant

/******************************************************************************/

			/*** 3. Cointegrating VAR Analysis of Real Money Balance: Johansens ML Approach ***/ 
  
	/* 3.1 VAR Order Selection */

varsoc mp y r, maxlag(4)        /* AIC: P=3 */

var mp y r, lags(1/3)
varstable

varlmar, mlag(3)
cap drop r_mp r_y r_r
predict r_mp, equation(#1) residuals 
predict r_y, equation(#2) residuals 
predict r_r, equation(#3) residuals 
corrgram r_mp, lags(32)
corrgram r_y, lags(32)
corrgram r_r, lags(32)

	/* 3.2 Johansens Likelihood Ratio Tests for the Rank of Cointegration */
	
vecrank mp y r, lags(3) trend(constant) max ic levela

	/* 3.3 ML Estimation of Vector Error Correction Model */
	
vec mp y r, lags(3) trend(constant) rank(1)

	/* 3.4 Granger Causality Test */
	
vec mp y r, lags(3) trend(constant) rank(1)
test[D_mp]: ld.y l2d.y ld.r l2d.r l._ce1 
display invchi2(5, 1-0.05)

	/* 3.5 Impulse Response Analysis and Forecast Error Variance Decomposition */
	
vec mp y r, lags(3) trend(constant) rank(1)
mat omega=e(omega)
mat list omega
scalar sd_mp=sqrt(omega[1,1])
scalar sd_y=sqrt(omega[2,2])
scalar sd_r=sqrt(omega[3,3])
scalar cor_mp_y=omega[1,2]/(sd_mp*sd_y)
scalar cor_mp_r=omega[1,3]/(sd_mp*sd_r)
scalar cor_y_r=omega[2,3]/(sd_y*sd_r)
scalar list sd_mp sd_y sd_r cor_mp_y cor_mp_r cor_y_r

/* Generalized IRF and FEVD */

/* shock to y */
irf set gir_gvd, replace
vec y mp r, lags(3) trend(constant) rank(1)
irf create order1, step(48)
preserve
use gir_gvd.irf, clear
replace oirf=oirf/sd_y if impulse=="y"
save gir_gvd.irf, replace
restore
irf table oirf, set(gir_gvd) irf(order1) impulse(y) response(mp y r)
irf graph oirf, set(gir_gvd) irf(order1) impulse(y) response(mp)
irf table fevd, set(gir_gvd) irf(order1) impulse(y) response(mp)	
irf graph fevd, set(gir_gvd) irf(order1) impulse(y) response(mp)


/* shock to r */
irf set gir_gvd, replace
vec r mp y, lags(3) trend(constant) rank(1)
irf create order2, step(48)
preserve
use gir_gvd.irf, clear
replace oirf=oirf/sd_r if impulse=="r"
save gir_gvd.irf, replace
restore
irf table oirf, set(gir_gvd) irf(order2) impulse(r) response(mp y r)
irf graph oirf, set(gir_gvd) irf(order2) impulse(r) response(mp)
irf table fevd, set(gir_gvd) irf(order2) impulse(r) response(mp)	
irf graph fevd, set(gir_gvd) irf(order2) impulse(r) response(mp)

/* shock to mp */
irf set gir_gvd, replace
vec mp y r, lags(3) trend(constant) rank(1)
irf create order3, step(48)
preserve
use gir_gvd.irf, clear
replace oirf=oirf/sd_mp if impulse=="mp"
save gir_gvd.irf, replace
restore
irf table oirf, set(gir_gvd) irf(order3) impulse(mp) response(mp y r)
irf graph oirf, set(gir_gvd) irf(order3) impulse(mp) response(mp) 
irf table fevd, set(gir_gvd) irf(order3) impulse(mp) response(mp)	
irf graph fevd, set(gir_gvd) irf(order3) impulse(mp) response(mp)


/******************************************************************************/

			/*** 4. Forecasting: VECM vs. VAR ***/ 


/* A. the VECM selected in part 3 */

/* 1-step-ahead */
vec mp y r if tin(,1982q4), lags(3) trend(constant) rank(1)
fcast compute h1_, step(1) replace
list date h1_mp if date=="1983Q1"

cap drop mp1_hat_table
gen mp1_hat_table = .
label var mp1_hat_table "VECM h=1"

set more off
local h=1					
local i=144                  
local l=`i'+`h'
local k=168                  

while `i' <=`k'-`h' {

qui vec mp y r in 1/`i', lags(3) trend(constant) rank(1)	
fcast compute hat_, step(`h') replace nose		
		
local j=`i'+`h' 
replace mp1_hat_table = hat_mp in `j'/`j'
local i=`i'+1
}
list date mp1_hat_table in `l'/`k'

tsline mp mp1_hat_table in `l'/`k'

cap drop mp1_fit 
cap drop mp_error mse rmspe mppe mappe
gen mp1_fit=mp1_hat_table
gen mp_error=mp1_fit - mp
replace mp_error=. if time<q(1983q1)
egen mse=mean(mp_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(mp_error/mp)
egen mappe=mean(abs(mp_error/mp))
list rmspe mppe mappe in 1/1

cap drop VECM_h1
gen VECM_h1=mp1_hat_table
label var VECM_h1 "VECM h=1"

/* 4-step-ahead */
vec mp y r if tin(,1982q4), lags(3) trend(constant) rank(1)
fcast compute h4_, step(4) replace
list date h4_mp if date=="1983Q1"

cap drop mp4_hat_table
gen mp4_hat_table = .
label var mp4_hat_table "VECM h=4"

set more off
local h=4					
local i=144                  
local l=`i'+`h'
local k=168                  

while `i' <=`k'-`h' {

qui vec mp y r in 1/`i', lags(3) trend(constant) rank(1)	
fcast compute hat_, step(`h') replace nose		
		
local j=`i'+`h' 
replace mp4_hat_table = hat_mp in `j'/`j'
local i=`i'+1
}
list date mp4_hat_table in `l'/`k'

tsline mp mp4_hat_table in `l'/`k'

cap drop mp4_fit 
cap drop mp_error mse rmspe mppe mappe
gen mp4_fit=mp4_hat_table
gen mp_error=mp4_fit - mp
replace mp_error=. if time<q(1983q1)
egen mse=mean(mp_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(mp_error/mp)
egen mappe=mean(abs(mp_error/mp))
list rmspe mppe mappe in 1/1

cap drop VECM_h4
gen VECM_h4=mp4_hat_table
label var VECM_h4 "VECM h=4"


/* B. Unrestricted VAR of first difference */

varsoc d.mp d.y d.r   /* AIC: p=2*/
var d.mp d.y d.r, lags(1/2)

/* I have trouble of forecasting the level variable using Unrestricted VAR of first difference */


/* C. Unrestricted VAR of first difference */

varsoc mp y r /* AIC: P=3*/
var mp y r, lags(1/3)

/* 1-step-ahead */
cap drop h1_hat_table 
gen h1_hat_table=.
cap drop yyy mp_fit mp_error mse rmspe mppe mappe

set more off
local h=1
local i=144
local l=`i'+`h'
local k=168
while `i'<=`k' -`h' {
quietly var mp y r in 1/`i', lags(1/3) 
fcast compute hat_, step(`h') replace nose

local j=`i'+`h'
replace h1_hat_table=hat_mp in `j'/`j'
local i=`i'+1
}
tsline mp h1_hat_table in `l'/`k'

cap drop yyy mp_fit
gen yyy=d.mp+l.mp
gen mp_fit =h1_hat_table+l.mp
tsline mp yyy mp_fit in `l'/`k' 

list date h1_hat_table

cap drop yyy mp_fit mp_error mse rmspe mppe mappe
gen yyy=mp
gen mp_fit =h1_hat_table
gen mp_error = mp_fit - yyy
egen mse=mean(mp_error^2)
gen rmspe=sqrt(mse)
egen mppe=mean(mp_error/yyy)
egen mappe=mean(abs(mp_error/yyy))
list rmspe mppe mappe in 1/1

cap drop VAR_level_h1 
gen VAR_level_h1=h1_hat_table
label var VAR_level_h1 "Unrestricted VAR (level) h=1"
tsline mp VAR_level_h1 if tin(1983q1,)	  

/* 4-step-ahead */
cap drop h4_hat_table 
gen h4_hat_table=.

set more off
local h=4
local i=144
local l=`i'+`h'
local k=168
while `i'<=`k' -`h' {
quietly var mp y r in 1/`i', lags(1/3)
fcast compute hat_, step(`h') replace nose

local j=`i'+`h'
replace h4_hat_table=hat_mp in `j'/`j'
local i=`i'+1
}
tsline mp h4_hat_table in `l'/`k'

list date h4_hat_table

cap drop yyy mp_fit mp_error mse rmspe mppe mappe
gen yyy=mp
gen mp_fit =h4_hat_table
gen mp_error=mp_fit -yyy
egen mse=mean(mp_error^2)
gen rmspe=sqrt(mse)
egen mppe=mean(mp_error/yyy)
egen mappe=mean(abs(mp_error/yyy))
list rmspe mppe mappe in 1/1

cap drop VAR_level_h4
gen VAR_level_h4=h4_hat_table
label var VAR_level_h4 "Unrestricted VAR (level) h=4"
tsline mp VAR_level_h4 if tin(1983q1,)

/* C. Compare */
tsline mp VECM_h1 VAR_level_h1 if tin(1983q1,)
tsline mp VECM_h4 VAR_level_h4 if tin(1983q1,) 
  
