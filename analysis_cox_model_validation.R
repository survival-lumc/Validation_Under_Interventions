#------------------------------
#------------------------------
#Obtains estimates used in the calculations of discrimination, calibration, and Brier score using validation data
#------------------------------
#------------------------------

#--------
#Values of L measured at time 0 - used later

if(scenario != 3) #note for scenario 3 L.baseline.dat is generated in the file dat_sim_cox_scenario3 incl (measurement) error
  {
    L.baseline.dat=dat.long.val$L[dat.long.val$visit==1]
  }

#------------------------------------
#------------------------------------
#ESTIMATED risks under the two treatment strategies: "always treated" (risk1), "never treated" (risk0)
#------------------------------------
#------------------------------------

#----
#estimated risks at fixed time horizons (t.hor) for each person 

#time horizons
t.hor=1:5

#risks under NEVER TREATED
risk0_exp=sapply(t.hor,FUN=function(x){1-exp(-(cumhaz.fun(x)*exp(cox.msm$coefficients["L.baseline"]*L.baseline.dat)))})

#risks under ALWAYS TREATED
risk1_exp=sapply(t.hor,FUN=function(x){1-exp(-(
  cumhaz.fun(min(x,1))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=1)*(cumhaz.fun(min(x,2))-cumhaz.fun(1))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=2)*(cumhaz.fun(min(x,3))-cumhaz.fun(2))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["Alag2"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=3)*(cumhaz.fun(min(x,4))-cumhaz.fun(3))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["Alag2"]+cox.msm$coefficients["Alag3"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)+
    (x>=4)*(cumhaz.fun(min(x,5))-cumhaz.fun(4))*exp(cox.msm$coefficients["A"]+cox.msm$coefficients["Alag1"]+cox.msm$coefficients["Alag2"]+cox.msm$coefficients["Alag3"]+cox.msm$coefficients["Alag4"]+cox.msm$coefficients["L.baseline"]*L.baseline.dat)
))})

#------------------------------
#------------------------------
#OBSERVERD risks under the two treatment strategies 'ALWAYS TREATED' and 'NEVER TREATED'
#------------------------------
#------------------------------

#-----------------
#'Observed' risks up to times 1:5, under the 'NEVER TREATED' strategy
#Obtained using censoring and weighting
#-----------------

#fit weights model for people untreated at time 0
wt.mod=glm(A~L,family = "binomial",data=dat.long.val[dat.long.val$Alag1==0,])

#predicted probability that A[t]=0 conditional on A[t-1]=0 and conditional on L[t]
pred.wt0=predict(wt.mod,type = "response",newdata = dat.long.val)
dat.long.val$wt0=1-pred.wt0

#Obtain the IPW at each time using cumulative product of 1/wt up to that time
dat.long.val = dat.long.val %>% group_by(id) %>% mutate(ipw0=1/cumprod(wt0))

#Now impose the artificial censoring when people deviate from the 'never treated' strategy

dat.long.val = dat.long.val %>%group_by(id) %>%mutate(A.baseline=first(A))
dat.long.val$in.dat.0 = (dat.long.val$A==0)

#---
#weighted Kaplan-Meier - using unstabilized weights
km.0=survfit(Surv(time,time.stop,event)~1, data=dat.long.val %>% filter(in.dat.0==1), weights = ipw0)
step.risk0.obs=stepfun(km.0$time,c(1,km.0$surv))#step function giving survival probability at any time

#estimated 'observed' risk at times 1:5
risk0_obs=1-step.risk0.obs(1:5)

#-----------------
#OBSERVED RISKS
#'Observed' risks up to times 1:5, under the 'ALWAYS TREATED' strategy
#Note that under our data generating mechanism people always continue treatment after they start
#i.e. there are no transitions from A=1 to A=0
#so we don't need time-dependent weights under the 'always treated' strategy
#-----------------

#---
#create data set which will be for people in the A=1 group at baseline

#fit weights model for people treated at time 0
wt.mod.baseline=glm(A~L,family = "binomial",data=dat.long.val[dat.long.val$visit==1,])


#predicted probability that A[0]=1 conditional on A[-1]=0 (which is true for everyone) and conditional on L[0]
pred.wt1.baseline=predict(wt.mod.baseline,type = "response",newdata = dat.long.val[dat.long.val$visit==1,])

#Obtain the IPW at each time (which is the same at each time here)
dat.long.val$wt1 = 0
dat.long.val$wt1[dat.long.val$visit==1] = 1/pred.wt1.baseline
dat.long.val = dat.long.val %>% group_by(id) %>% mutate(ipw1=sum(wt1))

#Now impose the artificial censoring when people deviate from the 'always treated' strategy
#Note once a person starts treatment they always continue, in this simulation study
dat.long.val$in.dat.1 = (dat.long.val$A.baseline==1)

#---
#weighted Kaplan-Meier
km.1=survfit(Surv(time,time.stop,event)~1,data=dat.long.val %>% filter(in.dat.1==1),weights = ipw1)
step.risk1.obs=stepfun(km.1$time,c(1,km.1$surv))#step function giving survival probability at any time

#estimated 'observed' risk at times 1:5
risk1_obs=1-step.risk1.obs(1:5)

#------------------------------------
#------------------------------------
#RESTRICTED SUBSET APPROACH
#indicators of being untreated or treated up to time 1-5
#------------------------------------
#------------------------------------

#indicator of being untreated (up to a given time point)
#Note that A is set to 0 after a person has the event (or is censored, though there is no censoring here)
subset0.1<-ifelse(dat.val$A.0==0,1,0)
subset0.2<-ifelse(dat.val$A.0==0 & dat.val$A.1==0,1,0)
subset0.3<-ifelse(dat.val$A.0==0 & dat.val$A.1==0 & dat.val$A.2==0,1,0)
subset0.4<-ifelse(dat.val$A.0==0 & dat.val$A.1==0 & dat.val$A.2==0 & dat.val$A.3==0,1,0)
subset0.5<-ifelse(dat.val$A.0==0 & dat.val$A.1==0 & dat.val$A.2==0 & dat.val$A.3==0 & dat.val$A.4==0,1,0)

subset0.matrix<-as.matrix(cbind(subset0.1,subset0.2,subset0.3,subset0.4,subset0.5))

#indicator of being treated (up to a given time point)
#Takes into account that A is set to 0 after a person has the event (or is censored, though there is no censoring here)
subset1.1<-ifelse(dat.val$A.0==1,1,0)
subset1.2<-ifelse((dat.val$A.0==1 & dat.val$A.1==1 & dat.val$T.obs>1)|(subset1.1==1 & dat.val$T.obs<=1),1,0)
subset1.3<-ifelse((dat.val$A.0==1 & dat.val$A.1==1 & dat.val$A.2==1 & dat.val$T.obs>2)|(subset1.2==1 & dat.val$T.obs<=2),1,0)
subset1.4<-ifelse((dat.val$A.0==1 & dat.val$A.1==1 & dat.val$A.2==1 & dat.val$A.3==1 & dat.val$T.obs>3)|(subset1.3==1 & dat.val$T.obs<=3),1,0)
subset1.5<-ifelse((dat.val$A.0==1 & dat.val$A.1==1 & dat.val$A.2==1 & dat.val$A.3==1 & dat.val$A.4==1 & dat.val$T.obs>4)|(subset1.4==1 & dat.val$T.obs<=4),1,0)

subset1.matrix<-as.matrix(cbind(subset1.1,subset1.2,subset1.3,subset1.4,subset1.5))

#-----------------
#RESTRICTED SUBSET APPROACH: OBSERVED RISKS
#Observed risks up to times 1:5, under the 'NEVER TREATED' strategy in the subset of people who are untreated up to that time
#-----------------

#---
#observed risks in the subset of people untreated up to each time

km.0.subset.time1=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset0.1==1,])
km.0.subset.time1.step<-stepfun(km.0.subset.time1$time,c(1,km.0.subset.time1$surv))

km.0.subset.time2=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset0.2==1,])
km.0.subset.time2.step<-stepfun(km.0.subset.time2$time,c(1,km.0.subset.time2$surv))

km.0.subset.time3=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset0.3==1,])
km.0.subset.time3.step<-stepfun(km.0.subset.time3$time,c(1,km.0.subset.time3$surv))

km.0.subset.time4=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset0.4==1,])
km.0.subset.time4.step<-stepfun(km.0.subset.time4$time,c(1,km.0.subset.time4$surv))

km.0.subset.time5=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset0.5==1,])
km.0.subset.time5.step<-stepfun(km.0.subset.time5$time,c(1,km.0.subset.time5$surv))

risk0_obs_subset=c(1-km.0.subset.time1.step(1),1-km.0.subset.time2.step(2),1-km.0.subset.time3.step(3),1-km.0.subset.time4.step(4),1-km.0.subset.time5.step(5))

#---
#observed risks in the subset of people treated up to each time

km.1.subset.time1=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset1.1==1,])
km.1.subset.time1.step<-stepfun(km.1.subset.time1$time,c(1,km.1.subset.time1$surv))

km.1.subset.time2=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset1.2==1,])
km.1.subset.time2.step<-stepfun(km.1.subset.time2$time,c(1,km.1.subset.time2$surv))

km.1.subset.time3=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset1.3==1,])
km.1.subset.time3.step<-stepfun(km.1.subset.time3$time,c(1,km.1.subset.time3$surv))

km.1.subset.time4=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset1.4==1,])
km.1.subset.time4.step<-stepfun(km.1.subset.time4$time,c(1,km.1.subset.time4$surv))

km.1.subset.time5=survfit(Surv(T.obs,D.obs)~1,data=dat.val[subset1.5==1,])
km.1.subset.time5.step<-stepfun(km.1.subset.time5$time,c(1,km.1.subset.time5$surv))

risk1_obs_subset=c(1-km.1.subset.time1.step(1),1-km.1.subset.time2.step(2),1-km.1.subset.time3.step(3),1-km.1.subset.time4.step(4),1-km.1.subset.time5.step(5))


