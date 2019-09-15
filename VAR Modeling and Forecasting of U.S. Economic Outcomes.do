/**********   VAR Modeling and Forecasting of U.S. Economic Outcomes **********/
/*************                      Chao Shi                      *************/

/******************************************************************************/

clear				/* drop all data and value labels from memory */
clear matrix
set memory 100m		/* reset the amount of memory allocated to Stata */ 

			/*** 0. Import Data ***/


/* A. Import data in .csv format */
insheet using "C:\Users\chris\Desktop\ECON825\US.csv", names comma clear
describe		/* check numeric vs. string variables */

/* B. Declare the data to be time series & specify the time variable: two alternatives */   			
gen time=quarterly(date,"yq")	/* quarterly() convert a string to a numeric time representation, 1960q1=0 */
tsset time, quarterly			/* gaps allowed in time series */
describe		/* double check if the numeric time variable is correct */
list time gdp85 gdpdef m1 govt              /* check for missing values */	


/* C. Log Transformation */
gen y=log(gdp85)
gen p=log(gdpdef)
gen m=log(m1)
gen x=log(govt)

codebook y p m x

/* D. Time Series Plot */

tsline y
tsline p 
tsline m 
tsline x
tsline d.y 
tsline d2.p 
tsline d.m 
tsline d.x

/******************************************************************************/

			/*** 1. VAR Modeling of (d:y; d2:p; d:m; d:x) ***/

/* 1.1 Standard VAR to (d:y; d2:p; d:m; d:x) */

/* A. Order Selection */
varsoc d.y d2.p d.m d.x, maxlag(8)		/* AIC: P=4, SBC: P=1 */
		
/* B. Estimation */
var d.y d2.p d.m d.x, lag(1)         /* P=1 */
var d.y d2.p d.m d.x, lag(1/4)         /* P=4 */

/* C. Diagnostics */
varlmar, mlag(4)		/* LM test for residual serial correlation */
varstable		/*  eigenvalues of the companion matrix, i.e., inverted roots of |phi(z)|=0 */             
varnorm			/* tests are applied to the orthogonalized residuals */

cap drop r_y r_p r_m r_x 
predict r_y, equation(#1) residuals
predict r_p, equation(#2) residuals
predict r_m, equation(#3) residuals
predict r_x, equation(#4) residuals

corrgram r_y, lags(32)	
corrgram r_p, lags(32)	
corrgram r_m, lags(32)	
corrgram r_x, lags(32)	


/* 1.2 Augmented VAR to (d:y; d2:p; d:m; d:x) */

/* Create Seasonal Dummies */
gen quarter=quarter(dofq(time))
gen sd_1=(quarter==1)
gen sd_2=(quarter==2)
gen sd_3=(quarter==3)
gen sd_4=(quarter==4)
list time quarter sd_* in 1/20, sep(4)


/* Create Structural Breaks Dummies (sbd_82q1 & sbd_87q1) */
cap drop sbd_82q1 sbd_87q1 
gen sbd_82q1=0
replace sbd_82q1=1 if time>=q(1982Q1)
gen sbd_87q1=0
replace sbd_87q1=1 if time>=q(1987Q1)
list time sbd_82q1 sbd_87q1

/* A. Order Selection */
varsoc d.y d2.p d.m d.x, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) maxlag(8)	/* AIC: P=3, SBC: P=0 */

/* B. Estimation  */
var d.y d2.p d.m d.x, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)

/* C. Diagnostics */

varlmar, mlag(4)	/* Pass */	
varstable   /*  eigenvalues of the companion matrix, i.e., inverted roots of |phi(z)|=0 */  
varnorm     /* tests are applied to the orthogonalized residuals */

cap drop r_y r_p r_m r_x 
predict r_y, equation(#1) residuals
predict r_p, equation(#2) residuals
predict r_m, equation(#3) residuals
predict r_x, equation(#4) residuals


corrgram r_y, lags(32)
corrgram r_p, lags(32)
corrgram r_m, lags(32)
corrgram r_x, lags(32)


/* 1.3 Testing the Effectiveness of Monetary and Fiscal Policies */
  
/* A. Granger Causality Test */
var d.y d2.p d.m d.x, exog(sd_1 sd_2 sd_3  sbd_82q1 sbd_87q1) lag(1/3)	/* augmented VAR */

/* if (d:y; d2:p) Granger-causes (d:m; d:x) */
test [D_m]ld.y [D_m]l2d.y [D_m]l3d.y [D_m]ld2.p [D_m]l2d2.p [D_m]l3d2.p [D_x]ld.y [D_x]l2d.y [D_x]l3d.y [D_x]ld2.p [D_x]l2d2.p [D_x]l3d2.p				  
display invchi2(12, 1-0.05)

/* if (d:m; d:x) Granger-causes (d:y; d2:p) */
test [D_y]ld.m [D_y]l2d.m [D_y]l3d.m [D_y]ld.x [D_y]l2d.x [D_y]l3d.x [D2_p]ld.m [D2_p]l2d.m [D2_p]l3d.m [D2_p]ld.x [D2_p]l2d.x [D2_p]l3d.x				  
display invchi2(12, 1-0.05)	

/* B. Test for Contemporaneous Effect via the Error Covariance Matrix */

var d.y d2.p d.m d.x, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1 ) lag(1/3)
scalar ll_sys=e(ll) /* system log-likelihood of the unrestricted model */

/* single equation log-likelihood of the restricted model */
reg d.y ld.y l2d.y l3d.y ld2.p l2d2.p l3d2.p ld.m l2d.m l3d.m ld.x l2d.x l3d.x sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1
scalar ll_dy=e(ll)
	
reg d2.p ld.y l2d.y l3d.y ld2.p l2d2.p l3d2.p ld.m l2d.m l3d.m ld.x l2d.x l3d.x sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1
scalar ll_d2p=e(ll)	

reg d.m ld.y l2d.y l3d.y ld2.p l2d2.p l3d2.p ld.m l2d.m l3d.m ld.x l2d.x l3d.x sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1
scalar ll_dm=e(ll)	

reg d.x ld.y l2d.y l3d.y ld2.p l2d2.p l3d2.p ld.m l2d.m l3d.m ld.x l2d.x l3d.x sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1
scalar ll_dx=e(ll)	  
			
scalar lrt=2*(ll_sys-ll_dy-ll_d2p-ll_dm-ll_dx)
scalar list ll_sys ll_dy ll_d2p ll_dm ll_dx lrt

display invchi2(6, 1-0.05)

/* 1.4 Impulse Response Analysis and Forecast Error Variance Decomposition */

cap drop dy d2p dm dx
gen dy=d.y
gen d2p=d2.p
gen dm=d.m
gen dx=d.x

 /* A: Orthogonalized IRF & FEVD */	
irf set oir_ovd, replace

var dy d2p dm dx, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)	/* augmented VAR */
						 
irf create order1, step(8) replace 

/* Orthogonalized IR of d.y and d2.p to a one-time unit shock to equation d.m */
irf table oirf, irf(order1) impulse(dm) response(dy d2p)						 
irf graph oirf, irf(order1) impulse(dm) response(dy d2p)

/* Orthogonalized VD of d.y and d2.p to a one-time unit shock to equation d.m */
irf table fevd, irf(order1) impulse(dm) response(dy d2p)						 
irf graph fevd, irf(order1) impulse(dm) response(dy d2p)


/* Orthogonalized IR of d.y and d2.p to a one-time unit shock to equation d.x */
irf table oirf, irf(order1) impulse(dx) response(dy d2p)						 
irf graph oirf, irf(order1) impulse(dx) response(dy d2p)

/* Orthogonalized VD of d.y and d2.p to a one-time unit shock to equation d.x */
irf table fevd, irf(order1) impulse(dx) response(dy d2p)						 
irf graph fevd, irf(order1) impulse(dx) response(dy d2p)


/* Check whether OVD over different innovations add up to 1 for forecast horizon h = 4 */

irf table fevd, irf(order1) impulse(dy d2p dm dx) response(dy)						 
irf graph fevd, irf(order1) impulse(dy d2p dm dx) response(dy)	

irf table fevd, irf(order1) impulse(dy d2p dm dx) response(d2p)						 
irf graph fevd, irf(order1) impulse(dy d2p dm dx) response(d2p)	

/* B. Generalized IRF & FEVD */

var dy d2p dm dx, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3) /* unrestricted model */	
mat omega=e(Sigma) /* error covariance matrix */
mat list omega
scalar sd_dy=sqrt(omega[1,1])   /* standard deviation of GDP innovations */
scalar sd_d2p=sqrt(omega[2,2])	/* standard deviation of inflation innovations */					 
scalar sd_dm=sqrt(omega[3,3])	/* standard deviation of monetary policy innovations */					 
scalar sd_dx=sqrt(omega[4,4])	/* standard deviation of fiscal policy innovations */					 

scalar list sd_dy sd_d2p sd_dm sd_dx

irf set gir_gvd, replace

					 
var dm dy d2p dx, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3) /* dm as first variable */	
irf create order1, step(8) replace

/* Generalized IR of d.y and d2.p to a one-time unit shock to equation d.m */								 
irf table oirf, irf(order1) impulse(dm) response(dy d2p)						 
irf graph oirf, irf(order1) impulse(dm) response(dy d2p)


/* Generalized VD of d.y and d2.p to a one-time unit shock to equation d.m */								 						 
irf table fevd, irf(order1) impulse(dm) response(dy d2p)						 
irf graph fevd, irf(order1) impulse(dm) response(dy d2p)						 
	

var dx dy d2p dm, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)	/* dx as first variable*/	
irf create order2, step(8) replace 
 
 /* Generalized IR of d.y and d2.p to a one-time unit shock to equation d.x */								 						 
irf table oirf, irf(order2) impulse(dx) response(dy d2p)						 
irf graph oirf, irf(order2) impulse(dx) response(dy d2p)
 
/* Generalized VD of d.y and d2.p to a one-time unit shock to equation d.x */								 						 
irf table fevd, irf(order2) impulse(dx) response(dy d2p)						 
irf graph fevd, irf(order2) impulse(dx) response(dy d2p)		



/* Check whether GVD over different innovations add up greater than 1 for forecast horizon h=4 */
				 
var dy dx d2p dm, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)	/* dy as first variable*/	
irf create order3, step(8) replace 
 
/* Generalized VD of d.y and d2.p to a one-time unit shock to equation d.y */								 						 
irf table fevd, irf(order3) impulse(dy) response(dy d2p)						 
irf graph fevd, irf(order3) impulse(dy) response(dy d2p)	


var d2p dx dy dm, exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)	/* d2p as first variable*/	
irf create order4, step(8) replace 

/* Generalized VD of d.y and d2.p to a one-time unit shock to equation d2.p */								 						 
irf table fevd, irf(order4) impulse(d2p) response(dy d2p)						 
irf graph fevd, irf(order4) impulse(d2p) response(dy d2p)


/******************************************************************************/

			/*** 2. Forecasting ***/

/* A. Dynammic forecasts */

var dy d2p dm dx if tin(,1986q4),exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3) 
fcast compute dyn_, step(20) replace				/* forecasts for 1987Q1-1991Q4 */
label var dyn_dy "dy, VARX(3), dynamic forecasts"
label var dyn_d2p "d2p, VARX(3), dynamic forecasts"

cap drop yyy y_fit y_error mse rmspe mppe mappe
gen yyy=dy
gen y_fit=dyn_dy
gen y_error=y_fit-yyy
replace y_error=. if time<tq(1987q1)
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1
tsline dy dyn_dy if tin(1987q1,)


cap drop yyy y_fit y_error mse rmspe mppe mappe
gen yyy=d2p
gen y_fit=dyn_d2p
gen y_error=y_fit-yyy
replace y_error=. if time<tq(1987q1)
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1
tsline d2p dyn_d2p if tin(1987q1,)


/* B. 1-step-ahead forecasts */

var dy d2p dm dx if tin(,1986q4),exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)
fcast compute h1_, step(1) replace		/* forecasts start in 1987Q1 */
list date h1_dy h1_d2p h1_dm h1_dx if date=="1987Q1"

cap drop y1_hat_table y2_hat_table 
gen y1_hat_table = .
gen y2_hat_table = .



set more off
local h=1					/* set the h for h-step-ahead-forcast */ 
local i=108					/* set last obs of estimation sample*/
local l=`i'+`h'
local k=128					/* set last obs of forecast sample */
while `i' <=`k'-`h' {

quietly var dy d2p dm dx in 1/`i',exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3) 
fcast compute hat_, step(`h') replace nose		
		
local j=`i'+`h' 

replace y1_hat_table = hat_dy in `j'/`j'
replace y2_hat_table = hat_d2p in `j'/`j'

local i=`i'+1
}
tsline dy y1_hat_table  in `l'/`k'
tsline d2p y2_hat_table in `l'/`k'



list date y1_hat_table y2_hat_table

cap drop yyy y_fit y_error mse rmspe mppe mappe
gen yyy=dy
gen y_fit=y1_hat_table
gen y_error=y_fit-yyy
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1


cap drop yyy y_fit y_error mse rmspe mppe mappe
gen yyy=d2p
gen y_fit=y2_hat_table
gen y_error=y_fit-yyy
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1


cap drop h1_dy h1_d2p h1_dm h1_dx
gen h1_dy=y1_hat_table
gen h1_d2p=y2_hat_table

label var h1_dy "dy, VARX(3), h=1"
label var h1_d2p "d2p, VARX(3), h=1"

tsline dy h1_dy if tin(1987q1,)
tsline d2p h1_d2p if tin(1987q1,)


/* C. 4-step-ahead forecasts */

var dy d2p dm dx if tin(,1986q4),exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3)
fcast compute h4_, step(4) replace		/* forecasts start in 1989Q4 */
list date h4_dy h4_d2p h4_dm h4_dx if date=="1987Q4"

cap drop y1_hat_table y2_hat_table 
gen y1_hat_table = .
gen y2_hat_table = .


set more off
local h=4					/* set the h for h-step-ahead-forcast */ 
local i=108					/* set last obs of estimation sample*/
local l=`i'+`h'
local k=128					/* set last obs of forecast sample */
while `i' <=`k'-`h' {

quietly var dy d2p dm dx in 1/`i',exog(sd_1 sd_2 sd_3 sbd_82q1 sbd_87q1) lag(1/3) 
fcast compute hat_, step(`h') replace nose		
		
local j=`i'+`h' 

replace y1_hat_table = hat_dy in `j'/`j'
replace y2_hat_table = hat_d2p in `j'/`j'

local i=`i'+1
}
tsline dy y1_hat_table in `l'/`k'
tsline d2p y2_hat_table in `l'/`k'

list date y1_hat_table y2_hat_table

cap drop yyy y_fit y_error mse rmspe mppe mappe
gen yyy=dy
gen y_fit=y1_hat_table
gen y_error=y_fit-yyy
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1

cap drop yyy y_fit y_error mse rmspe mppe mappe
gen yyy=d2p
gen y_fit=y2_hat_table
gen y_error=y_fit-yyy
egen mse=mean(y_error^2)
gen rmspe=sqrt(mse)      		
egen mppe=mean(y_error/yyy)
egen mappe=mean(abs(y_error/yyy))
list rmspe mppe mappe in 1/1

cap drop h4_dy h4_d2p
gen h4_dy=y1_hat_table
gen h4_d2p=y2_hat_table

label var h4_dy "dy, VARX(3), h=4"
label var h4_d2p "d2p,VARX(3), h=4"


tsline dy h4_dy if tin(1987q4,)
tsline d2p h4_d2p if tin(1987q4,)

/* D. forecast comparison */
tsline dy dyn_dy h1_dy h4_dy if tin(1987q1,)
tsline d2p dyn_d2p h1_d2p h4_d2p if tin(1987q1,)


