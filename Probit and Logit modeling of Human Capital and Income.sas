/***    Probit and Logit modeling of Human Capital and Income     ***/
/*********************** Chao Shi         ***************************/
/********************************************************************/

/*** import data ***/

libname project "C:\Users\Chao\Desktop\ECON804\project";

data data;
	set project.shi_chao;
run;

/*** create educ categorical variable ***/

data data1;
	set data;
	educ=.;
	if hgc=0 then educ=0;
	if 0<hgc<7 then educ=1;
	if hgc=7 then educ=2;
	if hgc>7 then educ=3;
run;

proc means data=data1 maxdec=2; var work age educ lnnonlabinc numchild; 
	title "TABLE 1: Descriptive Statistics";
run;

/******************** General Requirements *************************/

/*** 1. Randomly split the sample into two subsamples, with roughly 2/3 
of the observations in sample A and  the remaining 1/3 of the observations 
in sample B ***/

data data2;
	set data1;
	seed=100;
	random_order=ranuni(seed);
	if random_order>0.33333333 then sample=0;
	else sample=1;
run;

proc freq data=data2; table sample; 
	title "TABLE 2: Sample Split";
run;

data project.data2;
	set data2;
run;

/*** 4. Check the contingency table of the response variable against each 
categorical regressor to see if there are small/sparse cells 
(i.e., cell size less than 10) ***/

proc freq data=data2;
	tables work*educ/ nocol norow nopercent;
	title "TABLE 3: Contigency Table";
run;

/*** 5. Check if the regressors suffer multicollinearity (MC). 
Using sample A, report and interpret the Variance Inflation Factor, 
Condition Index, and Proportion of Variation ***/

data data_a;
	set data2;
	if sample=0;
run;

proc reg data=data_a plots=none;
	model work = age educ lnnonlabinc numchild  / vif tol collin collinoint;
	title "Table 4: MC Analysis 1";
run;

proc reg data=data_a plots=none;
	model work = age educ lnnonlabinc  / vif tol collin collinoint;
	title "Table 5: MC Analysis 2";
run;

/*** 6. Test the association between the response variable 
and each of the remaining categorical regressors ***/

proc freq data=data2;
	tables work*educ / chisq;
	title "Table 6: Association Test";
run;

/*** 7. Check the convergence status of the logit/probit model, 
make sure the convergence criterion is satisfied ***/

proc logistic data=data2 descending;
	model work = age educ lnnonlabinc;
	title "Convergence Check: Logit";
run;

proc logistic data=data2 descending;
	model work = age educ lnnonlabinc/ link=probit;
	title "Convergence Check: Probit";
run;


/******************** Part B.i. Binomial Logit Model *************************/

/*** 1. Estimate a binomial logit model. Report and interpret the estimation results ***/
/*** 2. Conduct the Hosmer-Lemeshow goodness-of-fit test ***/

proc logistic data=data_a descending;
	model work = age educ lnnonlabinc / lackfit;
	title "Table 7: Binomial logit model";
run;

/*** 3. Test the overall significance of each multi-level categorical regressor ***/

proc logistic data=data_a descending;
	class educ /param=effect;
	model work = age educ lnnonlabinc;
	contrast 'Joint Significance'  educ  1  0  0,
                                    educ  0  1  0,
                                    educ  0  0  1/estimate;
	title "Table 8: Multi-level categorical test";
run;

/*** 4. Odds Ratios: Compute and interpret the odds ratios for the following changes in the predictors ***/

/*** a. There is a change in one predictor at a time. For continuous variables,
consider a two-standard-error increase if the coefficient is positive, 
or a two-standard-error decrease if the coefficient is negative. 
For numerical but discretely-measured variables, such as years of school or number of children in a household, 
consider a one-unit increase if the coefficient is positive, or a one-unit decrease if the coefficient is negative. 
For categorical variables, odds ratios should be calculated by comparing each pair of outcome levels. 
Report the 95% confidence limits of these odds ratios. ***/

proc logistic data=data_a descending;
	model work = age educ lnnonlabinc /clodds=wald;
	units age=1 lnnonlabinc=2*sd;
	title "Table 9: Odds Ratio";
run;

data data_a1;
	set data_a;
	educ_0=0;
	educ_1=0;
	educ_2=0;
	educ_3=0;
	if educ=0 then educ_0=1;
	if educ=1 then educ_1=1;
	if educ=2 then educ_2=1;
	if educ=3 then educ_3=1;
run;

proc logistic data=data_a1 descending;
	model work = age educ_0 educ_1 educ_2 educ_3 lnnonlabinc /noint covb;
	title "Table 9: Odds Ratio";
run;

proc logistic data=data_a1 descending;
	model work = age educ_0 educ_1 educ_2 educ_3 lnnonlabinc /noint covb;
	units age=1 lnnonlabinc=2*sd;
	title "Table 9: Odds Ratio";
run;

/*** 5. Regression diagnosis: Examine the graphs of DIFDEV (called Deviance Deletion Difference in SAS output),
DIFCHISQ (called Chi-square Deletion Difference in SAS output), LEVERAGE, C and CBAR statistics ***/

proc logistic data=data_a descending; 
	model work = age educ lnnonlabinc / iplots;
	title "Table 10: Regression diagnosis";
run;

/*** 6. Evaluate the predictive power of the logit model 
Based on the parameter estimates reported in Part B.i.1, 
predict the probability of the “event” as well as the value of the response variable for each observation in sample B. 
The response variable is predicted to be an “event” if the predicted probability of event is greater 
than the predicted probability of non-event. These are the out-of-sample predictions. Print out the predicted ***/

data data_b;
	set data2;
	if sample=1;
run;

proc logistic data=data_a descending outmodel=logit_est noprint;
	model work = age educ lnnonlabinc;
run;

proc logistic inmodel=logit_est;
	score data=data_b out=logit_pred_data;
run;

data logit_pred_new;
	set logit_pred_data;
	id=_N_;
run;

proc print data=logit_pred_new (obs=10);
var age educ lnnonlabinc i_work P_1;title "Table 11: Predictive Power";run;


/******************** Part B.ii. Binomial Probit Model *************************/

/*** 1. Estimate a binomial probit model. Report and interpret the estimation results ***/
/*** 2. Conduct the Hosmer-Lemeshow goodness-of-fit test ***/

proc logistic data=data_a descending;
	model work = age educ lnnonlabinc / link=probit lackfit;
	output out=probit_out;
	title "Table 12: Binomial Probit Model";
run;

/*** 3. Odds: Compute the odds (not odds ratios) of the “event” for the following observations ***/
/*** a. Each categorical predictor is held at the level that will give the lowest probability of 
the “event”, holding all other predictors constant. Each numerical predictor takes a sequence of 
the following values: sample mean minus two standard deviations, sample mean, sample mean plus two 
standard deviations ***/

proc means data=probit_out;
	var age lnnonlabinc ;
	output out=std_means_out mean=avgage avglnnonlabinc  std=stdage stdlnnonlabinc;
	title "Table 13: Probit Odds";
run;

data mean_out_3a;
	set std_means_out;
	educ=0;
	do x=-2 to 2 by 2;
		age=avgage+(stdage*x);
		lnnonlabinc=avglnnonlabinc+(stdlnnonlabinc*x);
		output;
	end;
run;

proc logistic data=data_a outmodel=probit_est noprint;
	model work = age educ lnnonlabinc/link=probit;
run;

proc logistic inmodel=probit_est;
	score data=mean_out_3a out=probit_pred_3a;
run;

data work_odds_3a;
	set probit_pred_3a;
	odds=p_1/p_0;
run;

proc print data=work_odds_3a;
	var educ age lnnonlabinc odds;
	title "Table 13: Probit Odds";
run;

/*** b. Each categorical predictor is held at the level that will give the highest probability of 
the “event”, holding all other predictors constant. Each numerical predictor takes a sequence of 
the following values: sample mean minus two standard deviations, sample mean, sample mean plus two 
standard deviations ***/

data mean_out_3b;
	set std_means_out;
	educ=3;
	do x=-2 to 2 by 2;
		age=avgage+(stdage*x);
		lnnonlabinc=avglnnonlabinc+(stdlnnonlabinc*x);
		output;
	end;
run;

proc logistic inmodel=probit_est;
	score data=mean_out_3b out=probit_pred_3b;
run;

data work_odds_3b;
	set probit_pred_3b;
	odds=p_1/p_0;
run;

proc print data=work_odds_3b;
	var educ age lnnonlabinc odds;
	title "Table 13: Probit Odds";
run;

/*** 4 Odds Ratios: Based on the results in part 3, compute and interpret 
the odd ratios for the following simultaneous changes: each numerical predictor 
increases by the increment of two standard deviations, and each categorical predictor 
moves from the lowest-event-probability level to the highest-event-probability level ***/

/*** 4a ***/
data mean_out_4a;
	set std_means_out;
	do x=-2 to 0 by 2;
		age=avgage+(stdage*x);
		lnnonlabinc=avglnnonlabinc+(stdlnnonlabinc*x);
		output;
	end;
run;

data mean_out_4aa;
	set mean_out_4a;
	educ=.;
	if x=-2 then educ=0;
	if x=0 then educ=3;
run;

proc logistic inmodel=probit_est;
	score data=mean_out_4aa out=probit_pred_4a;
run;

data work_odds_4a;
	set probit_pred_4a;
	odds=p_1/p_0;
run;

proc print data=work_odds_4a;
	var educ age lnnonlabinc odds;
	title "Table 14: Probit Odds Ratio";
run;

/*** 4b ***/
data mean_out_4b;
	set std_means_out;
	do x=0 to 2 by 2;
		age=avgage+(stdage*x);
		lnnonlabinc=avglnnonlabinc+(stdlnnonlabinc*x);
		output;
	end;
run;

data mean_out_4bb;
	set mean_out_4b;
	educ=.;
	if x=0 then educ=0;
	if x=2 then educ=3;
run;

proc logistic inmodel=probit_est;
	score data=mean_out_4bb out=probit_pred_4b;
run;

data work_odds_4b;
	set probit_pred_4b;
	odds=p_1/p_0;
run;

proc print data=work_odds_4b;
	var educ age lnnonlabinc odds;
	title "Table 14: Probit Odds Ratio";
run;

/*** 5. Regression diagnosis: Examine the graphs of DIFDEV (called Deviance Deletion Difference in SAS output),
DIFCHISQ (called Chi-square Deletion Difference in SAS output), LEVERAGE, C and CBAR statistics ***/

proc logistic data=data_a descending; 
	model work = age educ lnnonlabinc / link=probit iplots;
	title "Table 15: Probit Regression diagnosis";
run;

/*** 6. Evaluate the predictive power of the logit model 
Based on the parameter estimates reported in Part B.ii.1, 
predict the probability of the “event” as well as the value of the response variable for each observation in sample B. 
The response variable is predicted to be an “event” if the predicted probability of event is greater 
than the predicted probability of non-event. These are the out-of-sample predictions. Print out the predicted ***/

proc logistic data=data_a descending outmodel=probit_est noprint;
	model work = age educ lnnonlabinc/link=probit;
run;

proc logistic inmodel=probit_est;
	score data=data_b out=probit_pred_data;
run;

data probit_pred_new;
	set probit_pred_data;
	id=_N_;
run;

proc print data=probit_pred_new (obs=10);
	var age educ lnnonlabinc i_work P_1;
	title "Table 16: Probit predictive power";
run;






