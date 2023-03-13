#------------------------------
#------------------------------
#An additive hazards MSM (using IPTW) is fitted to the counterfactual data, conditional on baseline L0.
#The MSM is correctly specified. No weights are required.
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
dat.true.long$L.baseline <- dat.true.long$L.0

#-----------------
#MSM-IPTW analysis using a correctly-specified additive hazards model
#-----------------

ah.true.p0=aalen(Surv(time,time.stop,event)~A+L.baseline,data=dat.true.long[dat.true.long$time==0,],n.sim=0)
ah.true.p1=aalen(Surv(time,time.stop,event)~A+Alag1+L.baseline,data=dat.true.long[dat.true.long$time==1,],n.sim=0)
ah.true.p2=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+L.baseline,data=dat.true.long[dat.true.long$time==2,],n.sim=0)
ah.true.p3=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+L.baseline,data=dat.true.long[dat.true.long$time==3,],n.sim=0)
ah.true.p4=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4+L.baseline,data=dat.true.long[dat.true.long$time==4,],n.sim=0)

maxrow.p0=dim(ah.true.p0$cum)[1]
maxrow.p1=dim(ah.true.p1$cum)[1]
maxrow.p2=dim(ah.true.p2$cum)[1]
maxrow.p3=dim(ah.true.p3$cum)[1]
maxrow.p4=dim(ah.true.p4$cum)[1]

ah.true.stepfun.int=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                       c(0,ah.true.p0$cum[,2],
                         ah.true.p0$cum[maxrow.p0,2]+ah.true.p1$cum[-1,2],
                         ah.true.p0$cum[maxrow.p0,2]+ah.true.p1$cum[maxrow.p1,2]+ah.true.p2$cum[-1,2],
                         ah.true.p0$cum[maxrow.p0,2]+ah.true.p1$cum[maxrow.p1,2]+ah.true.p2$cum[maxrow.p2,2]+ah.true.p3$cum[-1,2],
                         ah.true.p0$cum[maxrow.p0,2]+ah.true.p1$cum[maxrow.p1,2]+ah.true.p2$cum[maxrow.p2,2]+ah.true.p3$cum[maxrow.p3,2]+ah.true.p4$cum[-1,2]))

ah.true.stepfun.A=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                     c(0,ah.true.p0$cum[,3],
                       ah.true.p0$cum[maxrow.p0,3]+ah.true.p1$cum[-1,3],
                       ah.true.p0$cum[maxrow.p0,3]+ah.true.p1$cum[maxrow.p1,3]+ah.true.p2$cum[-1,3],
                       ah.true.p0$cum[maxrow.p0,3]+ah.true.p1$cum[maxrow.p1,3]+ah.true.p2$cum[maxrow.p2,3]+ah.true.p3$cum[-1,3],
                       ah.true.p0$cum[maxrow.p0,3]+ah.true.p1$cum[maxrow.p1,3]+ah.true.p2$cum[maxrow.p2,3]+ah.true.p3$cum[maxrow.p3,3]+ah.true.p4$cum[-1,3]))

ah.true.stepfun.Alag1=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.true.p0$cum)[1]),
                           ah.true.p1$cum[-1,4],
                           ah.true.p1$cum[maxrow.p1,4]+ah.true.p2$cum[-1,4],
                           ah.true.p1$cum[maxrow.p1,4]+ah.true.p2$cum[maxrow.p2,4]+ah.true.p3$cum[-1,4],
                           ah.true.p1$cum[maxrow.p1,4]+ah.true.p2$cum[maxrow.p2,4]+ah.true.p3$cum[maxrow.p3,4]+ah.true.p4$cum[-1,4]))

ah.true.stepfun.Alag2=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.true.p0$cum)[1]),
                           rep(0,dim(ah.true.p1$cum)[1]-1),
                           ah.true.p2$cum[-1,5],
                           ah.true.p2$cum[maxrow.p2,5]+ah.true.p3$cum[-1,5],
                           ah.true.p2$cum[maxrow.p2,5]+ah.true.p3$cum[maxrow.p3,5]+ah.true.p4$cum[-1,5]))

ah.true.stepfun.Alag3=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.true.p0$cum)[1]),
                           rep(0,dim(ah.true.p1$cum)[1]-1),
                           rep(0,dim(ah.true.p2$cum)[1]-1),
                           ah.true.p3$cum[-1,6],
                           ah.true.p3$cum[maxrow.p3,6]+ah.true.p4$cum[-1,6]))

ah.true.stepfun.Alag4=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.true.p0$cum)[1]),
                           rep(0,dim(ah.true.p1$cum)[1]-1),
                           rep(0,dim(ah.true.p2$cum)[1]-1),
                           rep(0,dim(ah.true.p3$cum)[1]-1),
                           ah.true.p4$cum[-1,7]))

ah.true.stepfun.L=stepfun(c(ah.true.p0$cum[,1], ah.true.p1$cum[-1,1],ah.true.p2$cum[-1,1],ah.true.p3$cum[-1,1],ah.true.p4$cum[-1,1]),
                     c(0,ah.true.p0$cum[,4],
                       ah.true.p0$cum[maxrow.p0,4]+ah.true.p1$cum[-1,5],
                       ah.true.p0$cum[maxrow.p0,4]+ah.true.p1$cum[maxrow.p1,5]+ah.true.p2$cum[-1,6],
                       ah.true.p0$cum[maxrow.p0,4]+ah.true.p1$cum[maxrow.p1,5]+ah.true.p2$cum[maxrow.p2,6]+ah.true.p3$cum[-1,7],
                       ah.true.p0$cum[maxrow.p0,4]+ah.true.p1$cum[maxrow.p1,5]+ah.true.p2$cum[maxrow.p2,6]+ah.true.p3$cum[maxrow.p3,7]+ah.true.p4$cum[-1,8]))

#-----------------
#Estimated 'true' risks under the treatment strategies 'ALWAYS TREATED' and 'NEVER TREATED'
#-----------------

#Values of L measured at time 0 - used later
L.baseline.dat=dat.cf$L.0

#----
#estimated true risks at fixed time horizons (t.hor) for each person 

#time horizons
t.hor=1:5

#risks under NEVER TREATED
risk0_true=sapply(FUN=function(h){1-exp(-(ah.true.stepfun.int(h)+ah.true.stepfun.L(h)*L.baseline.dat))},t.hor)

#risks under ALWAYS TREATED
risk1_true=sapply(FUN=function(h){1-exp(-(ah.true.stepfun.int(h)+ah.true.stepfun.A(h)+ah.true.stepfun.Alag1(h)+ah.true.stepfun.Alag2(h)+ah.true.stepfun.Alag3(h)+ah.true.stepfun.Alag4(h)+ah.true.stepfun.L(h)*L.baseline.dat))},t.hor)
