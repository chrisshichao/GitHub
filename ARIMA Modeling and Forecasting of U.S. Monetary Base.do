/********     ARIMA Modeling and Forecasting of U.S. Monetary Base    *********/
/*************                      Chao Shi                      *************/

/******************************************************************************/
/*** Import data from csv format ***/
insheet using "C:\Users\chris\Desktop\ECON825\US_IFS_2006.CSV", names comma clear
describe			
list time base

/*** Declare the data to be time series & specify the time variable: two alternatives ***/    
  /*Alternative 1*/
gen time_1=q(1957q1)+_n-1		
format time_1 %tq				
tsset time_1	

 /*Alternative 2*/
gen time_2=quarterly(time,"yq")	
tsset time_2, quarterly			

describe
list time time_1 time_2			

save "C:\Users\chris\Desktop\ECON825\US_IFS_2006.dta", replace
use "C:\Users\chris\Desktop\ECON825\US_IFS_2006.dta", clear 	

/*** Summary Statistics & Time Series Plots ***/
  /* summary statistics */
codebook base
  /* time series plot */
tsline base
/******************************************************************************/

/***** 1. Unit Root Tests *****/

/*** 1.1 Standard ADF Test ***/

/* 1. Plot data. */ 

gen yyy=log(base)
label var yyy "log(base)"
tsline yyy
tsline d.yyy
tsline d2.yyy
gen trend=_n

/* 
   2. Specify the ADF regression equation. 
   3. Select the proper lag order (order of augmentation) of the ADF regression
      according to AIC.
   4. Report the ADF t-statistics and the 5% critical values.
   5. Decide whether or not to reject H0.   
*/

/* Schwert (1989), pmax=int[12*{(T+1)/100}^0.25]=14 */
scalar T=191  /*9 missing value(8 at the beginning, 1 at end)*/
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=8 /* 8 of missing observations at the beginning of the sample */

/* level */
scalar diff=0
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg
  
/*with intercept and trend in the regression*/
reg yyy trend l.yyy l(1/14)d.yyy
estat ic
matrix y_order=r(S)
reg yyy trend l.yyy l(1/13)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/12)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/11)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/10)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/9)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/8)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/7)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/6)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/5)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/4)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/3)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy l(1/2)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy ld.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))

matlist y_order		/* AIC: -1094.655 df=16 P=16-3=13 */
dfuller yyy, trend lag(13) /*not reject at 5% significance level*/


/* first-difference */
scalar diff=1
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg
 
/*from the plot, no trend for the first difference (may have contant)*/
reg d.yyy ld.yyy l(1/14)d2.yyy
estat ic
matrix dy_order=r(S)
reg d.yyy ld.yyy l(1/13)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/12)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/11)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/10)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/9)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/8)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/7)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/6)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/5)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/4)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/3)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy l(1/2)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy ld2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))

matlist dy_order		/* AIC:  -1090.143 df=16 P=16-2=14 */
dfuller d.yyy, lag(14)  /*reject at 5% significance level*/

/* 6. Draw conclusions on the order of integration of yyy. */
/* AIC: yyy ~ I(1) with drift */	


/*** 1.2 ADF Test with Seasonal Dummies ***/

gen quarter=quarter(dofq(time_2))
gen sd_1=(quarter==1)
gen sd_2=(quarter==2)
gen sd_3=(quarter==3)
gen sd_4=(quarter==4)

list time_2 quarter sd_* in 1/20, sep(4)

/* 
   2. Specify the ADF regression equation. 
   3. Select the proper lag order (order of augmentation) of the ADF regression
      according to AIC.
   4. Report the ADF t-statistics and the 5% critical values.
   5. Decide whether or not to reject H0.   
*/

scalar T=191
scalar pmax=int(12*((T+1)/100)^0.25)
scalar miss=8

/* level */
scalar diff=0
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/14)d.yyy
estat ic
matrix y_order=r(S)
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/13)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/12)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/11)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/10)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/9)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/8)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/7)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/6)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/5)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/4)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/3)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/2)d.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 ld.yyy if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))
reg yyy trend l.yyy sd_1 sd_2 sd_3 if _n>p_reg
estat ic
matrix y_order=(y_order \ r(S))

matlist y_order		/* AIC: -1105.646 df=14 P=14-6=8 */

reg yyy trend l.yyy sd_1 sd_2 sd_3 l(1/8)d.yyy  /* not reject at 5% significance level*/


/* first-difference */
scalar diff=1
scalar p_reg=pmax+1+diff+miss
scalar list pmax diff p_reg

reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/14)d2.yyy
estat ic
matrix dy_order=r(S)
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/13)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/12)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/11)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/10)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/9)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/8)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/7)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/6)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/5)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/4)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/3)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 l(1/2)d2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 ld2.yyy if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))
reg d.yyy ld.yyy sd_1 sd_2 sd_3 if _n>p_reg
estat ic
matrix dy_order=(dy_order \ r(S))

matlist dy_order		/* AIC: -1099.044 df=12 P=12-5=7 */
reg d.yyy sd_1 sd_2 sd_3 l.d.yyy l(1/7).d2.yyy /*(no trend) reject at 5% significance level*/

/* 6. Draw conclusions on the order of integration of yyy. */
/* AIC: yyy ~ I(1) with drift */	

/******************************************************************************/

/***** 2. Seasonal Unit Root Test *****/

/*If the absolute value of test statistic is greater than the critical value, 
declare statistical significance and reject the null hypothesis.*/

hegy4 yyy, lag(8) det(strend) /*pi2,3,4, are significant              (include constant, trend and seasonal dummies)*/	
hegy4 d.yyy, lag(7) det(seas) /*Pi1,2,3, are significant              (include constant, trend and seasonal dummies)*/


/******************************************************************************/

/***** 3. ARIMA Modeling *****/

/* Both ADF tests with or w/o seasonal dummies indicating yyy is I(1)with drift*/
/* MA-AC AR-PAC */
/* p<0.05 - do not pass Q */

/*** 3.1. ARMA Modeling: d.yyy ***/ 
tsline d.yyy
corrgram d.yyy, lags(24)		
ac d.yyy, lags(24)	level(95)	
pac d.yyy, lags(24) level(95)	

/*1. Standard ARMA model */	
arima yyy, arima(4,1,4) /*coeffs. not significant;not pass unit root test*/
arima yyy, arima(3,1,3) /*not pass Q test and unit root test*/
/*
arima yyy, arima(4,1,2) /*not pass unit root test*/
arima yyy, arima(4,1,1) /*not pass unit root test*/
*/
arima yyy, arima(2,1,2) /*coeffs. not significant;not pass Q test and unit root test*/
arima yyy, arima(1,1,1) /*not pass Q test and unit root test*/

/*diagnostic checking*/
cap drop res
predict res, residuals
tsline res
corrgram res, lags(24)
wntestq res, lags(24)		/* Q-tests for white noise */
armaroots		/* check the stationarity & invertibilty condition, common AR & MA roots */


/* 2. Additive Seasonal ARMA model */

arima d.yyy, ar(1 4) ma(1 4) /*not pass unit root test*/

/*diagnostic checking*/
cap drop res
predict res, residuals
tsline res
corrgram res, lags(24)
wntestq res, lags(24)		/* Q-tests for white noise */
armaroots		/* check the stationarity & invertibilty condition, common AR & MA roots */


/* 3. Pure SARMA model */

arima d.yyy, mar(1,4) mma(1,4) /*not pass unit root test*/


/*diagnostic checking*/
cap drop res
predict res, residuals
tsline res
corrgram res, lags(24)
wntestq res, lags(24)		/* Q-tests for white noise */
armaroots		/* check the stationarity & invertibilty condition, common AR & MA roots */



/* 4. Multiplicative seasonal ARIMA models */

arima yyy, arima(1,1,1) sarima(1,0,1,4) /*not pass unit root test*/


   /*diagnostic checking*/
cap drop res
predict res, residuals
tsline res
corrgram res, lags(24)
wntestq res, lags(24)		/* Q-tests for white noise */
armaroots		/* check the stationarity & invertibilty condition, common AR & MA roots */


/*** 3.2  Structureal ARMAX models ***/

/*create residual variable*/
reg d.yyy sd_1 sd_2 sd_3
cap drop res
predict res, r

tsline res
corrgram res, lags(24)		
ac res, lags(40) level(95)	/*clear cutoff at 2 or 4*/
pac res, lags(40) level(95) /*clear cutoff at 2 or 4*/

 
arima d.yyy sd_1 sd_2 sd_3, ar(1/4)              /*coefficients for the 1st and 3rd lags are not significant*/
arima d.yyy sd_1 sd_2 sd_3, ar(2 4)              /* keep*/ 
arima d.yyy sd_1 sd_2 sd_3, ma(1/4)              /*coefficients for the 1st and 3rd lags are not significant*/
arima d.yyy sd_1 sd_2 sd_3, ma(2 4)              /* keep*/
arima yyy sd_1 sd_2 sd_3, arima(4,1,4)           /*some roots near unit circle and some coefficients insignificant*/
arima yyy sd_1 sd_2 sd_3, arima(2,1,4)           /*keep*/
arima yyy sd_1 sd_2 sd_3, arima(4,1,2)           /*some coefficients insignificant*/


/*diagnostic checking*/
cap drop res
predict res, residuals
tsline res
corrgram res, lags(24)
wntestq res, lags(24)		
armaroots

/*** model selection: use common sample ***/
quietly arima d.yyy sd_1 sd_2 sd_3, ar(2 4)
estat ic
quietly arima d.yyy sd_1 sd_2 sd_3, ma(2 4) 
estat ic
quietly arima yyy sd_1 sd_2 sd_3, arima(2,1,4)           /*Final Model: smallest AIC*/
estat ic

/******************************************************************************/

/***** 4. Forecasting *****/

/*** 4.1 In-sample tted values (estimation sample 1957Q1-2006Q2) ***/
     
arima yyy sd_1 sd_2 sd_3, arima(2,1,4)        /*full sample*/

cap drop y_fit
cap drop y_fit dy_fit	
predict y_fit, y 			
label var y_fit "In-sample fitted values of levels"
list time y_fit
tsline yyy y_fit if tin(1999q1,)

cap drop y_error mse rmspe mppe mappe
gen y_error=y_fit-yyy
replace y_error=. if time_2<tq(1999q1)		/* compute the forecast evaluation statistics for post-1999Q1 sample */
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1


/*** 4.2 Dynamic forecasts (xed estimation sample 1957Q1-1998Q4, varying forecast horizon) ***/
/* estimation sample is fixed as (1957q1, 1998q4), forecast horizon ranges from 1 to 30 */
/* dynamic forecasts start with 1999q1 */

arima yyy sd_1 sd_2 sd_3, arima(2,1,4) 

cap drop y_dyn 	
predict y_dyn, y dynamic(tq(1999q1))		/* forecasts start in 1999Q1 */			
label var y_dyn "dynamic forecasts of levels"
list time y_dyn if time_2>=tq(1999q1)
tsline yyy y_dyn if tin(1999q1,)

cap drop y_error mse rmspe mppe mappe
gen y_error=y_dyn-yyy
replace y_error=. if time_2<tq(1999q1)
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1

/*** 4.3 One-step-ahead forecasts (expanding estimation sample that begins with 1957Q1-1998Q4, forecast horizon xed at h=1) ***/
/* estimation sample starts with (1957q1, 1998q4), then expands with one additonal observation at a time, forecast horizon is fixed as h=1 */
/* reestimate the model with each additional observation */

arima yyy sd_1 sd_2 sd_3, arima(2,1,4) 
cap drop y_h1
predict y_h1, y 
label var y_h1 "1-step-ahead forecasts of levels"
list time y_h1 if time_2==tq(1999q1)		/* forecasts start in 1999Q1 */


cap drop y_hat
gen y_hat = .
cap drop y_hat_table
gen y_hat_table = .
label var y_hat_table "yyy, ARIMA, h=1"

set more off
local h = 1                  /* set the h for h-step-ahead-forcast */ 
local i=168                  /* set last obs of estimation sample */
local l=`i'+`h'
local k=198                  /* set last obs of forecast sample */

while `i' <=`k'-`h' {


quietly arima yyy sd_1 sd_2 sd_3 in 1/`i', arima(2,1,4) 


cap drop y_hat	
predict y_hat, y 
		
local j=`i'+`h' 
*	list time y_hat in `j'/`j'
replace y_hat_table = y_hat in `j'/`j'

local i=`i'+1
}

list time yyy y_hat_table in `l'/`k'
tsline yyy y_hat_table in `l'/`k'

cap drop y_error mse rmspe mppe mappe
gen y_error=y_hat_table-yyy
egen mse= mean(y_error^2)
gen rmspe= sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1

cap drop y_h1
gen y_h1=y_hat_table


/*** 4. Four-step-ahead forecasts (expanding estimation sample that begins with 1957Q1-1998Q4, forecast horizon xed at h=4) ***/
/* estimation sample starts with (1957q1, 1998q4), then expands with one additonal observation at a time, forecast horizon is fixed as h=4 */
/* reestimate the model with each additional observation */

arima yyy sd_1 sd_2 sd_3, arima(2,1,4) 

cap drop y_h4
predict y_h4, y dynamic(tq(1999q1))
label var y_h4 "4-step-ahead forecasts of levels"
list time y_h4 if time_2==tq(1999q4)		/* forecasts start in 1999Q4 */


cap drop y_hat
gen y_hat = .
cap drop y_hat_table
gen y_hat_table = .
label var y_hat_table "yyy, ARIMA, h=4"

set more off
local h = 4                  /* set the h for h-step-ahead-forcast */ 
local i=168                  /* set last obs of estimation sample */
local t0=tq(1998q4)			/* set ending date of estimation sample */
local l=`i'+`h'
local k=199                  /* set last obs of forecast sample */

while `i' <=`k'-`h' {

quietly arima yyy sd_1 sd_2 sd_3 in 1/`i', arima(2,1,4) 
cap drop y_hat	
predict y_hat, y dynamic(`t0'+1)
		
local j=`i'+`h' 
*	list time y_hat in `j'/`j'
replace y_hat_table = y_hat in `j'/`j'

local i=`i'+1
local t0=`t0'+1
}
list time yyy y_hat_table in `l'/`k'
tsline yyy y_hat_table in `l'/`k'

cap drop y_error mse rmspe mppe mappe
gen y_error=y_hat_table-yyy
egen mse= mean(y_error^2)
gen rmspe= sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1

cap drop y_h4
gen y_h4=y_hat_table


/* put all four forecast lines togethwer with the real data */
tsline yyy y_fit y_dyn y_h1 y_h4 if tin(1999q1,2006q3)
