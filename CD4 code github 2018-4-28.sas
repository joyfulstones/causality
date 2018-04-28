libname ddiddc "c:\data";

data ddiddc;
set ddiddc.ddiddc;
stratum=stratum-1;
gender=gender-1;
hemobl=hemobl-12;
id=seq;
trt=randgrp-1;
run;
* Change the data structure: each CD4 count as one row;
data one;
set ddiddc;
array cd4_all {11} cd4bl cd402 cd404 cd406 cd408 cd410 cd412 cd414 cd416 cd418 cd420;
aa=1;
stoptime=t2death/30;
status=death;
do i=1 to 11;
	cd4=cd4_all{i};
	month=2*(i-1);
	output;
end;
run;
data one2;
set one;
if cd4 ne .;
run;

data two;
set one2;
by id;
if first.id;
aa=2;
bb=1;
run;

* Get all the death event times;
data three;
set two;
if status=1;
run;
* Calculate the quantiles of the recurrent event time;
proc univariate data=three noprint;
var stoptime; 
output out=quant pctlpts=0 10 20 30 40 50 60 70 80 90 100 pctlpre=q; 
run;
data quant;
set quant;
bb=1;
run;
* Merge data with the quantiles;
data four;
merge two quant;
by bb;
run;

* Calculate the duration in each quantile interval, together with the indicator of event in each interval;
data five;
set four;
array quant {11} (0 2 4 6 8 10 12 14 16 18 20);
array dur {10} dur1-dur10;
array median {10} median1-median10;

array event {10} event1-event10;

do i=1 to 10;
	dur{i}=0;
	event{i}=0;
	median{i}=0;
end;

do i=2 to 11;
	if stoptime<=quant{i} then do;
		dur{i-1}=stoptime-quant{i-1};
		event{i-1}=status;
		median{i-1}=quant{i-1}+dur{i-1}/2; /* Get the median of each interval */
		lastmedian=median{i-1};
		i=11;
	end;
	else do;
		dur{i-1}=quant{i}-quant{i-1};
		median{i-1}=quant{i-1}+dur{i-1}/2; /* Get the median of each interval */
		lastmedian=median{i-1};
	end;
end;
run;
data six;
set one2 five;
log_cd4=log(cd4+1);
run;
data seven;
set six;
if log_cd4=. then log_cd4=stoptime;
year=month/12;
run;
proc sort data=seven;
by id aa;
run;
* Model II;
proc nlmixed data=seven qpoints=5;
parms h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
beta0=3.6 beta1=-.13 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23
alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
gamma1=-.4 vara=1.75 vare=.3;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara >=0;

base_haz=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
cum_base_haz=h1 * dur1 + h2 * dur2 + h3 * dur3 + h4 * dur4 + h5 * dur5 + h6 * dur6 + 
h7 * dur7 + h8* dur8 +h9 * dur9 + h10 * dur10;

if aa=1 then do;
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + a ;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;
if aa=2 then do;
	mu2= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl + gamma1 * a ;
	loglik2=-exp(mu2) * cum_base_haz;
	if status=0 then loglik=loglik2;
	if status=1 then loglik= log(base_haz) + mu2 + loglik2; 	/*log likelihood for failure */
end;
model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
predict (log_cd4-mu1)*(aa=1) out=ddiddc.res2; /* Output the residual in the mixed model */
predict (aa=2) *((status=1)+loglik2) out=ddiddc.mres2;	/* Output the martingale residual */
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;
data eight;
set seven;
array cdfour {10} cdfour1-cdfour10;
cdfour1=log(cd4bl+1);
cdfour2=log(cd402+1);
cdfour3=log(cd404+1);
cdfour4=log(cd406+1);
cdfour5=log(cd408+1);
cdfour6=log(cd410+1);
cdfour7=log(cd412+1);
cdfour8=log(cd414+1);
cdfour9=log(cd416+1);
cdfour10=log(cd418+1);
* Last value carried forward imputation;
* The baseline CD4 is always measured;
do i=2 to 10;
	if cdfour{i}=. then cdfour{i}=cdfour{i-1};
	last_cd4=cdfour{i};
end;
run;
* Model I by Proc NLMIXED;

proc nlmixed data=eight qpoints=5 maxiter=2000;

parms 	h1=0.02 h2=0.04 h3=0.04 h4=0.04 h5=0.05 h6=0.1 h7=0.06 h8=0.07 h9=0.05 h10=0.1
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23
		alpha1=-.3 alpha2=0 alpha3=.7 alpha4=0 alpha5=-.3 eta=-.3;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;

if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + a;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k});
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl  ;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl + eta * cd4_current ;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;
model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
run;

data eightb;
set eight;
if aa=1;
run;
data nine;
set eightb;
by id;
retain last_t last_logcd4;
if first.id then do;
	last_t=month;
	last_logcd4=log(cd4+1);
end;
else do;
		t1=last_t;
		t2=month;
		last_t=t2;
		logcd4=last_logcd4;
		last_logcd4=log(cd4+1);
		event=0;
		output;
end;
if last.id then do;
		t1=last_t;
		t2=stoptime;
		logcd4=last_logcd4;
		event=status;
		output;
	
end;
run;
* Model I by Proc PHRED, used in Table 1;

proc phreg data=nine;
   model (T1,T2)*event(0)=trt gender prevoi stratum hemobl logcd4;
run;	
* Model I-B by Proc PHRED, used in Table 1;

proc phreg data=nine;
   model (T1,T2)*event(0)=trt gender prevoi stratum hemobl logcd4 logcd4*trt;
run;	
* Model III;

proc nlmixed data=eight qpoints=5;

parms 	h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
		alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23
		gamma1=-.1 eta=0
		vara=1.75  vare=.3 ;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + a;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k});
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl +gamma1 * a ;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl + eta * cd4_current + gamma1 * a ;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;

model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 
predict (log_cd4-mu1)*(aa=1) out=ddiddc.cd4res3; /* Output the residual in the mixed model */
predict (aa=2) *((status=1)+loglik2) out=ddiddc.cd4mres3;	/* Output the martingale residual */

run;

* Model III with squared random effect;

proc nlmixed data=eight qpoints=5;

parms 	h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
		alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23
		gamma1=-.1 gamma2=0 eta=0
		vara=1.75  vare=.3 ;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + a;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k});
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl +gamma1 * a +gamma2 * a*a;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl + eta * cd4_current + gamma1 * a +gamma2 * a*a;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;

model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 

run;

* Model IV, with interaction of  a_i and Y(t);
* Used in the paper;

proc nlmixed data=eight qpoints=5;

parms 	h1=0.02 h2=0.02 h3=0.02 h4=0.02 h5=0.02 h6=0.02 h7=0.02 h8=0.02 h9=0.02 h10=0.02
		alpha1=-.28 alpha2=.3 alpha3=1.5 alpha4=-.16 alpha5=-.5
		beta0=3.6 beta1=-.1 beta2=-.8 beta3=0 beta4=-1.3 beta5=.12 beta6=.23
		gamma1=-.1 eta=0 delta=0
		vara=1.75  vare=.3 ;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;
array cdfour {10} cdfour1-cdfour10;


if aa=1 then do;	/* likelihood for CD4 measurements */
	mu1= beta0 +  beta1 * trt + beta2 * year + beta3* gender + beta4*prevoi + beta5*stratum + beta6 * hemobl + a;
	loglik=-.5*(log_cd4-mu1)**2/vare-.5*log(vare);
end;

if aa=2 then do;
	
		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent CD4 measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * cdfour{k} +  delta* cdfour{k} * a);
		end;

		mu3= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl +gamma1 * a ;

		loglik2=-exp(mu3) * sum2;

		if status=0 then loglik=loglik2;		/* for censoring event */
		if status=1 then do; 				/* for death event */
			base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
				h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
			cd4_current=cdfour1 * event1 + cdfour2 * event2 + cdfour3 * event3 + cdfour4 * event4 + cdfour5 * event5 + 
				cdfour6 * event6 + cdfour7 * event7 + cdfour8* event8 +cdfour9 * event9 + cdfour10 * event10;
			mu4= alpha1 * trt + alpha2 * gender + alpha3* prevoi + alpha4* stratum + alpha5 * hemobl 
				+ eta * cd4_current +   delta* cd4_current *a + gamma1 * a ;	
			loglik=loglik2 + mu4 + log(base_haz_d);
		end;
end;

model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 

predict (log_cd4-mu1)*(aa=1) out=ddiddc.cd4res4; /* Output the residual in the mixed model */
predict (aa=2) *((status=1)+loglik2) out=ddiddc.cd4mres4;	/* Output the martingale residual */

run;
data _null_;
set ddiddc.cd4res4;
file "c:\data\cd4res4.txt";
put id pred trt gender prevoi stratum year hemobl;
where pred ne 0;
run;
data _null_;
set ddiddc.cd4mres4;
file "c:\data\cd4mres4.txt";
put id pred trt gender prevoi stratum hemobl;
where pred ne 0;
run;
