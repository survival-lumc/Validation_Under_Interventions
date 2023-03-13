#------------------------------
#------------------------------
#A Cox MSM (using IPTW) is fitted to the development data, conditional on baseline L0.
#------------------------------
#------------------------------

#------------------
#Estimate stabilized weights
#The weights take into account that once a person starts treatment they always continue (in the simulated data)
#------------------

#Denominator of weights
wt.mod=glm(A~L,family="binomial",data=dat.long.dev[dat.long.dev$Alag1==0,])
pred.wt=predict(wt.mod,type = "response",newdata = dat.long.dev)
dat.long.dev$wt=ifelse(dat.long.dev$A==1,pred.wt,1-pred.wt)
dat.long.dev$wt=ifelse(dat.long.dev$Alag1==1,1,dat.long.dev$wt)
dat.long.dev$wt.cum=ave(dat.long.dev$wt,dat.long.dev$id,FUN=cumprod)

#Numerator of stabilised weights
wt.mod.num.L=glm(A~L.baseline*as.factor(visit),family="binomial",data=dat.long.dev[dat.long.dev$Alag1==0,])
pred.wt.num.L=predict(wt.mod.num.L,type = "response",newdata = dat.long.dev)
dat.long.dev$wt.num.L=ifelse(dat.long.dev$A==1,pred.wt.num.L,1-pred.wt.num.L)
dat.long.dev$wt.num.L=ifelse(dat.long.dev$Alag1==1,1,dat.long.dev$wt.num.L)
dat.long.dev$wt.cum.num.L=ave(dat.long.dev$wt.num.L,dat.long.dev$id,FUN=cumprod)

#Stabilized weights
dat.long.dev$ipw.s.L=dat.long.dev$wt.cum.num.L/dat.long.dev$wt.cum

#-----------------
#MSM-IPTW analysis using stabilized weights
#-----------------

cox.msm=coxph(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4+L.baseline,data=dat.long.dev,weights = dat.long.dev$ipw.s.L)

#baseline cumulative hazard
cumhaz=basehaz(cox.msm,centered=F)$hazard
event.times=basehaz(cox.msm,centered=F)$time

#step function giving baseline cumulative hazard at any time
cumhaz.fun=stepfun(event.times,c(0,cumhaz))
