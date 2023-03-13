#------------------------------
#------------------------------
#An additive hazards MSM (using IPTW) is fitted to the development data, conditional on baseline L0.
#Note that the MSM is correctly specified
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
#MSM-IPTW analysis using a correctly-specified additive hazards model, using stabilized weights
#-----------------

ah.p0=aalen(Surv(time,time.stop,event)~A+L.baseline,data=dat.long.dev[dat.long.dev$time==0,],n.sim=0,weights = dat.long.dev[dat.long.dev$time==0,]$ipw.s.L)
ah.p1=aalen(Surv(time,time.stop,event)~A+Alag1+L.baseline,data=dat.long.dev[dat.long.dev$time==1,],n.sim=0,weights = dat.long.dev[dat.long.dev$time==1,]$ipw.s.L)
ah.p2=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+L.baseline,data=dat.long.dev[dat.long.dev$time==2,],n.sim=0,weights = dat.long.dev[dat.long.dev$time==2,]$ipw.s.L)
ah.p3=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+L.baseline,data=dat.long.dev[dat.long.dev$time==3,],n.sim=0,weights = dat.long.dev[dat.long.dev$time==3,]$ipw.s.L)
ah.p4=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4+L.baseline,data=dat.long.dev[dat.long.dev$time==4,],n.sim=0,weights = dat.long.dev[dat.long.dev$time==4,]$ipw.s.L)

maxrow.p0=dim(ah.p0$cum)[1]
maxrow.p1=dim(ah.p1$cum)[1]
maxrow.p2=dim(ah.p2$cum)[1]
maxrow.p3=dim(ah.p3$cum)[1]
maxrow.p4=dim(ah.p4$cum)[1]

ah.stepfun.int=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                       c(0,ah.p0$cum[,2],
                         ah.p0$cum[maxrow.p0,2]+ah.p1$cum[-1,2],
                         ah.p0$cum[maxrow.p0,2]+ah.p1$cum[maxrow.p1,2]+ah.p2$cum[-1,2],
                         ah.p0$cum[maxrow.p0,2]+ah.p1$cum[maxrow.p1,2]+ah.p2$cum[maxrow.p2,2]+ah.p3$cum[-1,2],
                         ah.p0$cum[maxrow.p0,2]+ah.p1$cum[maxrow.p1,2]+ah.p2$cum[maxrow.p2,2]+ah.p3$cum[maxrow.p3,2]+ah.p4$cum[-1,2]))

ah.stepfun.A=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                     c(0,ah.p0$cum[,3],
                       ah.p0$cum[maxrow.p0,3]+ah.p1$cum[-1,3],
                       ah.p0$cum[maxrow.p0,3]+ah.p1$cum[maxrow.p1,3]+ah.p2$cum[-1,3],
                       ah.p0$cum[maxrow.p0,3]+ah.p1$cum[maxrow.p1,3]+ah.p2$cum[maxrow.p2,3]+ah.p3$cum[-1,3],
                       ah.p0$cum[maxrow.p0,3]+ah.p1$cum[maxrow.p1,3]+ah.p2$cum[maxrow.p2,3]+ah.p3$cum[maxrow.p3,3]+ah.p4$cum[-1,3]))

ah.stepfun.Alag1=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.p0$cum)[1]),
                           ah.p1$cum[-1,4],
                           ah.p1$cum[maxrow.p1,4]+ah.p2$cum[-1,4],
                           ah.p1$cum[maxrow.p1,4]+ah.p2$cum[maxrow.p2,4]+ah.p3$cum[-1,4],
                           ah.p1$cum[maxrow.p1,4]+ah.p2$cum[maxrow.p2,4]+ah.p3$cum[maxrow.p3,4]+ah.p4$cum[-1,4]))

ah.stepfun.Alag2=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.p0$cum)[1]),
                           rep(0,dim(ah.p1$cum)[1]-1),
                           ah.p2$cum[-1,5],
                           ah.p2$cum[maxrow.p2,5]+ah.p3$cum[-1,5],
                           ah.p2$cum[maxrow.p2,5]+ah.p3$cum[maxrow.p3,5]+ah.p4$cum[-1,5]))

ah.stepfun.Alag3=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.p0$cum)[1]),
                           rep(0,dim(ah.p1$cum)[1]-1),
                           rep(0,dim(ah.p2$cum)[1]-1),
                           ah.p3$cum[-1,6],
                           ah.p3$cum[maxrow.p3,6]+ah.p4$cum[-1,6]))

ah.stepfun.Alag4=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                         c(0,rep(0,dim(ah.p0$cum)[1]),
                           rep(0,dim(ah.p1$cum)[1]-1),
                           rep(0,dim(ah.p2$cum)[1]-1),
                           rep(0,dim(ah.p3$cum)[1]-1),
                           ah.p4$cum[-1,7]))

ah.stepfun.L=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                     c(0,ah.p0$cum[,4],
                       ah.p0$cum[maxrow.p0,4]+ah.p1$cum[-1,5],
                       ah.p0$cum[maxrow.p0,4]+ah.p1$cum[maxrow.p1,5]+ah.p2$cum[-1,6],
                       ah.p0$cum[maxrow.p0,4]+ah.p1$cum[maxrow.p1,5]+ah.p2$cum[maxrow.p2,6]+ah.p3$cum[-1,7],
                       ah.p0$cum[maxrow.p0,4]+ah.p1$cum[maxrow.p1,5]+ah.p2$cum[maxrow.p2,6]+ah.p3$cum[maxrow.p3,7]+ah.p4$cum[-1,8]))
