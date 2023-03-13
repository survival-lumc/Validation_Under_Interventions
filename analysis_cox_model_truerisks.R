#------------------------------
#------------------------------
#A Cox MSM is fitted to the counterfactual data, conditional on baseline L0
#No weights are required.
#This is used to calculate the true individual risks
#------------------------------
#------------------------------

#-----------------
#The counterfactual data has two components - never treated an always untreated 
#Combine these two  data sets into a combined data set 
#In this combined data set each individual has two copies
#-----------------

#---
#first put the counterfactual data into long frmat

dat.true.A1<-dat.cf[,c("id","L.0","D.A1","T.A1")]
dat.true.long.A1<-survSplit(formula=Surv(T.A1,D.A1) ~.,data=dat.true.A1,cut=c(0:5),start="time",end="time.stop",event="event")
dat.true.long.A1$A<-1
dat.true.long.A1$Alag1<-1*I(dat.true.long.A1$time>=1)
dat.true.long.A1$Alag2<-1*I(dat.true.long.A1$time>=2)
dat.true.long.A1$Alag3<-1*I(dat.true.long.A1$time>=3)
dat.true.long.A1$Alag4<-1*I(dat.true.long.A1$time>=4)

dat.true.A0<-dat.cf[,c("id","L.0","D.A0","T.A0")]
dat.true.long.A0<-survSplit(formula=Surv(T.A0,D.A0) ~.,data=dat.true.A0,cut=c(0:5),start="time",end="time.stop",event="event")
dat.true.long.A0$A<-0
dat.true.long.A0$Alag1<-0
dat.true.long.A0$Alag2<-0
dat.true.long.A0$Alag3<-0
dat.true.long.A0$Alag4<-0

#---
#Then combine into a single data set

dat.true.long<-rbind(dat.true.long.A1,dat.true.long.A0)

#-----------------
#MSM-IPTW analysis using a Cox model
#-----------------

cox.msm.true=coxph(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4+L.0,data=dat.true.long,control = coxph.control(timefix = FALSE))

#baseline cumulative hazard
cumhaz=basehaz(cox.msm.true,centered=F)$hazard
event.times=basehaz(cox.msm.true,centered=F)$time

#step function giving baseline cumulative hazard at any time
cumhaz.fun=stepfun(event.times,c(0,cumhaz))

#-----------------
#Estimated 'true' risks under the two treatment strategies 'ALWAYS TREATED' and 'NEVER TREATED'
#-----------------

#Values of L measured at time 0 - used later
L.baseline.dat=dat.cf$L.0

#----
#estimated true risks at fixed time horizons (t.hor) for each person 

#time horizons
t.hor=1:5

#risks under NEVER TREATED
risk0_true=sapply(t.hor,FUN=function(x){1-exp(-(cumhaz.fun(x)*exp(cox.msm.true$coefficients["L.0"]*L.baseline.dat)))})

#risks under ALWAYS TREATED
risk1_true=sapply(t.hor,FUN=function(x){1-exp(-(
  cumhaz.fun(min(x,1))*exp(cox.msm.true$coefficients["A"]+cox.msm.true$coefficients["L.0"]*L.baseline.dat)+
    (x>=1)*(cumhaz.fun(min(x,2))-cumhaz.fun(1))*exp(cox.msm.true$coefficients["A"]+cox.msm.true$coefficients["Alag1"]+cox.msm.true$coefficients["L.0"]*L.baseline.dat)+
    (x>=2)*(cumhaz.fun(min(x,3))-cumhaz.fun(2))*exp(cox.msm.true$coefficients["A"]+cox.msm.true$coefficients["Alag1"]+cox.msm.true$coefficients["Alag2"]+cox.msm.true$coefficients["L.0"]*L.baseline.dat)+
    (x>=3)*(cumhaz.fun(min(x,4))-cumhaz.fun(3))*exp(cox.msm.true$coefficients["A"]+cox.msm.true$coefficients["Alag1"]+cox.msm.true$coefficients["Alag2"]+cox.msm.true$coefficients["Alag3"]+cox.msm.true$coefficients["L.0"]*L.baseline.dat)+
    (x>=4)*(cumhaz.fun(min(x,5))-cumhaz.fun(4))*exp(cox.msm.true$coefficients["A"]+cox.msm.true$coefficients["Alag1"]+cox.msm.true$coefficients["Alag2"]+cox.msm.true$coefficients["Alag3"]+cox.msm.true$coefficients["Alag4"]+cox.msm.true$coefficients["L.0"]*L.baseline.dat)
))})
