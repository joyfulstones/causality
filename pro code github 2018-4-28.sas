PROC IMPORT OUT= WORK.pro 
            DATAFILE= "c:\data\pro.csv" 
            DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
data a;
set pro;
by id;
fup=stop;
if last.id;
keep id fup;
run;
data one;
merge pro a;
by id;
aa=1;
run;
data two;
set one;
by id;
if last.id;
aa=2;
bb=1;
run;

* Get all the death event times;
data three;
set two;
if death=1;
run;
* Calculate the quantiles of the death time;
proc univariate data=three noprint;
var fup; 
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
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array dur {10} dur1-dur10;

array event {10} event1-event10;

do i=1 to 10;
	dur{i}=0;
	event{i}=0;
end;

do i=2 to 11;
	if fup<=quant{i} then do;
		dur{i-1}=fup-quant{i-1};
		event{i-1}=death;
		i=11;
	end;
	else do;
		dur{i-1}=quant{i}-quant{i-1};
	end;
end;
run;
data five2;
set one five;
run;
proc sort data=five2;
by id aa;
run;

* Model II, used in the paper;
proc nlmixed data=five2 qpoints=5;
parms h1=0.2 h2=0.2 h3=0.2 h4=0.2 h5=0.2 h6=0.2 h7=0.2 h8=0.2 h9=0.2 h10=0.2
beta0=76 beta1=-6.5 beta2=1.7
alpha1=0.1 gamma1=0 vara=359 vare=344.76;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara >=0;

base_haz=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
cum_base_haz=h1 * dur1 + h2 * dur2 + h3 * dur3 + h4 * dur4 + h5 * dur5 + h6 * dur6 + 
h7 * dur7 + h8* dur8 +h9 * dur9 + h10 * dur10;

if aa=1 then do;
	mu1= beta0 +  beta1 * trt + beta2 * start + a ;
	loglik=-.5*(pro-mu1)**2/vare-.5*log(vare);
end;
if aa=2 then do;
	mu2= alpha1 * trt + gamma1 * a ;
	loglik2=-exp(mu2) * cum_base_haz;
	if death=0 then loglik=loglik2;
	if death=1 then loglik= log(base_haz) + mu2 + loglik2; 	/*log likelihood for failure */
end;
model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;

run;

proc univariate data=one;
var pro;
run;
data quant2;
set quant;
aa=1;
run;
data six;
merge one quant2;
by aa;
pro_c=pro-79;
run;
* For each Pro observational time, we need to calculate the duration in each quantile interval of death time;
data six2;
set six;
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array dur {10} dur1-dur10;

last_start=start;
do i=1 to 10;
	dur{i}=0;
end;

do i=2 to 11;
	if stop>quant{i} then do;
		if last_start<= quant{i} then do;
			dur{i-1}=quant{i}-last_start;
			last_start=quant{i};
		end;
	end;
	else do;
		dur{i-1}=stop - last_start;
		i=11;
	end;
end;
run;	

data six3;
set six2;
array quant {11} q0 q10 q20 q30 q40 q50 q60 q70 q80 q90 q100;
array event {10} event1-event10;
by id;
lastid=0;
if last.id then do;
	lastid=1;
	do i=1 to 10;
		event{i}=0;
	end;

	do i=2 to 11;
		if stop<=quant{i} then do;
			event{i-1}=death;
			i=11;
		end;
	end;
end;
run;
* Create status=0 for each repeated measures except the last.id;
data six4;
set six3;
by id;

if last.id then do;
		status=death;
		output;
end;
else do;
		status=0;
		output;
end;
run;

* Model I by Proc PHREG, used in Table 2;

proc phreg data=six4;
   model (start,stop)*status(0)=trt pro ;
run;	
ods graphics off;
proc phreg data=six4;
   model (start,stop)*status(0)=trt pro_c ;
run;	
* Model I-B by Proc PHREG;

proc phreg data=six4;
   model (start,stop)*status(0)=trt pro_c trt* pro_c ;
run;	
* Model 0;
proc mixed data=six3 covtest;
class id;
model pro = trt | start /s;
random id;
run;

* Model III, used in Table 2;

proc nlmixed data=six3 qpoints=5;
parms 	h1=3 h2=3 h3=3 h4=2 h5=1.5 h6=2 h7=2 h8=2 h9=2.5 h10=5
		beta0=76 beta1=-6.8 beta2=1.5 
		alpha1=-.14 gamma1=0 eta=0
		vara=359  vare=344;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;

/* likelihood for repeated measures */
	mu1= beta0 +  beta1 * trt + beta2 * start + a ;
	loglik1=-.5*(pro-mu1)**2/vare-.5*log(vare);

		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * pro);
		end;

  	    mu3= alpha1 * trt +gamma1 * a ;
		loglik2=-exp(mu3) * sum2;

		loglik=loglik1 + loglik2;

		if lastid=1 then do;

			if death=1 then do; 				/* for death event */
				base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
					h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
				mu4= alpha1 * trt + eta * pro + gamma1 *a ;	
				loglik=loglik + mu4 + log(base_haz_d);
			end;
		end;
	
model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 

run;
* Model III with squared random effect;

proc nlmixed data=six3 qpoints=5;
parms 	h1=3 h2=3 h3=3 h4=2 h5=1.5 h6=2 h7=2 h8=2 h9=2.5 h10=5
		beta0=76 beta1=-6.8 beta2=1.5 
		alpha1=-.14 gamma1=0 gamma2=0 eta=0
		vara=359  vare=344;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;

/* likelihood for repeated measures */
	mu1= beta0 +  beta1 * trt + beta2 * start + a;
	loglik1=-.5*(pro-mu1)**2/vare-.5*log(vare);

		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * pro);
		end;

  	    mu3= alpha1 * trt +gamma1 * a +gamma2 * a*a;
		loglik2=-exp(mu3) * sum2;

		loglik=loglik1 + loglik2;

		if lastid=1 then do;

			if death=1 then do; 				/* for death event */
				base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
					h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
				mu4= alpha1 * trt + eta * pro + gamma1 *a +gamma2 * a*a;	
				loglik=loglik + mu4 + log(base_haz_d);
			end;
		end;
	
model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
run;

* Model IV with interaction of random effect and Y(t);
* Used in the paper;

proc nlmixed data=six3 qpoints=5;
parms 	h1=3 h2=3 h3=3 h4=2 h5=1.5 h6=2 h7=2 h8=2 h9=2.5 h10=5
		beta0=76 beta1=-6.8 beta2=1.5 
		alpha1=-.14 gamma1=0 eta=0 delta=0 
		vara=359  vare=344;
bounds h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 vara  >=0;

array dur {10} dur1-dur10;
array baseh {10} h1-h10;

/* likelihood for repeated measures */
	mu1= beta0 +  beta1 * trt + beta2 * start + a;
	loglik1=-.5*(pro-mu1)**2/vare-.5*log(vare);

		sum2=0;
		do k=1 to 10;
			/* cumulative baseline hazard for time dependent measure */
			sum2=sum2 + baseh{k} * dur{k} * exp(eta * pro +   delta * pro * a);
		end;

  	    mu3= alpha1 * trt +gamma1 * a ;
		loglik2=-exp(mu3) * sum2;

		loglik=loglik1 + loglik2;

		if lastid=1 then do;

			if death=1 then do; 				/* for death event */
				base_haz_d=h1 * event1 + h2 * event2 + h3 * event3 + h4 * event4 + h5 * event5 + 
					h6 * event6 + h7 * event7 + h8* event8 +h9 * event9 + h10 * event10;
				mu4= alpha1 * trt + eta * pro  + delta * pro * a + gamma1 *a ;	
				loglik=loglik + mu4 + log(base_haz_d);
			end;
		end;
	
model id ~ general(loglik);
random a  ~ normal([0], [vara]) subject=id;
ods output ParameterEstimates=est1 FitStatistics=fit1; 

predict (pro-mu1) out=prores4; /* Output the residual in the mixed model */
predict (loglik2) out=promres4;	/* Output the martingale residual */
run;

data _null_;
set prores4;
file "c:\data\prores4.txt";
res=pred;
put id res trt start;
where pred ne 0;
run;
proc sql;
create table mres4 as
select distinct
id,
sum(pred) as mres
from promres4
group by id;
run;
/*
Proc summary data=promres4 nway;
    class id;
    var pred;
    output out = mres4 sum=;
run;
*/
* Calculate the Martingal residual: note due to the time varying covariate, we need to sum up loglik2 
in each row for the same subject;
data mres4;
merge two mres4;
by id;
res=death+pred;
run;
data _null_;
set mres4;
file "y:\liulei\paper\paper20\promres4.txt";
put id res trt time;
*where pred ne 0;
run;

