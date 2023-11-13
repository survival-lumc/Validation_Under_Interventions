#------------------------------------
#------------------------------------
#ESTIMATED mean risks at times 1:5 under the two treatment strategies: "always treated" (risk1), "never treated" (risk0)
#------------------------------------
#------------------------------------

calib_risk0[i,,1]=sapply(1:5,FUN=function(x){mean(risk0_exp[,x])})
calib_risk1[i,,1]=sapply(1:5,FUN=function(x){mean(risk1_exp[,x])})

#------------------------------------
#------------------------------------
#OBSERVED risks at times 1:5 under the two treatment strategies: "always treated" (risk1), "never treated" (risk0)
#------------------------------------
#------------------------------------

calib_risk0[i,,2]=risk0_obs
calib_risk1[i,,2]=risk1_obs

#------------------------------------
#------------------------------------
#TRUE risks at times 1:5 under the two treatment strategies: "always treated" (risk1), "never treated" (risk0)
#------------------------------------
#------------------------------------

#---
#never treated
km.true.0=survfit(Surv(T.A0,D.A0)~1,data=dat.cf)
km.true.0.step<-stepfun(km.true.0$time,c(1,km.true.0$surv))
calib_risk0[i,,3]=1-km.true.0.step(1:5)

#---
#always treated
km.true.1=survfit(Surv(T.A1,D.A1)~1,data=dat.cf)
km.true.1.step<-stepfun(km.true.1$time,c(1,km.true.1$surv))
calib_risk1[i,,3]=1-km.true.1.step(1:5)

#------------------------------------
#------------------------------------
#Mean predicted and observed risks within categories of predicted risks under the two treatment strategies: "always treated", "never treated"
#------------------------------------
#------------------------------------

#-----
#generate expected risk groups
#USING DATA SET-SPECIFIC CUT POINTS

cut_points_0=c(0,quantile(risk0_exp[,5],probs=seq(0.1,0.9,0.1)),1)
cut_points_1=c(0,quantile(risk1_exp[,5],probs=seq(0.1,0.9,0.1)),1)

risk0_group=cut(risk0_exp[,5],breaks=cut_points_0,include.lowest = T,labels = F)
dat.risk0.grp=data.frame(id=1:n,risk0_group)
dat.0.merge=merge(dat.risk0.grp, dat.long.val %>% filter(in.dat.0==1), by="id")

risk1_group=cut(risk1_exp[,5],breaks=cut_points_1,include.lowest = T,labels = F)
dat.risk1.grp=data.frame(id=1:n,risk1_group)
dat.1.merge=merge(dat.risk1.grp,dat.long.val %>% filter(in.dat.1==1),by="id")

#-----
#mean predicted risk within each decile

calib_risk0_group[i,,1]<-sapply(FUN=function(x){mean(risk0_exp[,5][risk0_group==x])},1:10)
calib_risk1_group[i,,1]<-sapply(FUN=function(x){mean(risk1_exp[,5][risk1_group==x])},1:10)

#-----
#observed risk in each decile of risk - "never treated"

km.0.grp=survfit(Surv(time,time.stop,event)~strata(risk0_group),data=dat.0.merge, weights = dat.0.merge$ipw0)

risk0_obs_grp<-rep(NA,10)
for(k in 1:10){
  #This line is because occasionally there are no individuals in a given group using the cut points. 
  group.exists<-paste0("strata(risk0_group)=risk0_group=",k)%in%summary(km.0.grp,cens=T)$strata
  
  if(group.exists){
    step.grp=stepfun(km.0.grp$time[summary(km.0.grp,cens=T)$strata==paste0("strata(risk0_group)=risk0_group=",k)],
                        c(1,km.0.grp$surv[summary(km.0.grp,cens=T)$strata==paste0("strata(risk0_group)=risk0_group=",k)]))
    risk0_obs_grp[k]<-1-step.grp(5)
  }
}
calib_risk0_group[i,,2]<-risk0_obs_grp

#-----
#observed risk in each decile of risk - "always treated"

km.1.grp=survfit(Surv(time,time.stop,event)~strata(risk1_group),data=dat.1.merge,weights = dat.1.merge$ipw1)

risk1_obs_grp<-rep(NA,10)
for(k in 1:10){
  #This line is because occasionally there are no individuals in a given group using the cut points. 
  group.exists<-paste0("strata(risk1_group)=risk1_group=",k)%in%summary(km.1.grp,cens=T)$strata
  
  if(group.exists){
    step.grp=stepfun(km.1.grp$time[summary(km.1.grp,cens=T)$strata==paste0("strata(risk1_group)=risk1_group=",k)],
                     c(1,km.1.grp$surv[summary(km.1.grp,cens=T)$strata==paste0("strata(risk1_group)=risk1_group=",k)]))
    risk1_obs_grp[k]<-1-step.grp(5)
  }
}
calib_risk1_group[i,,2]<-risk1_obs_grp

#------------------------------------
#------------------------------------
#TRUE risks within categories of predicted risks under the two treatment strategies: "always treated", "never treated"
#------------------------------------
#------------------------------------

#-----
#true risk in each decile of risk - "never treated"

km.0.grp=survfit(Surv(T.A0,D.A0)~strata(risk0_group),data=dat.cf)

risk0_true_grp<-rep(NA,10)
for(k in 1:10){
  #This line is because occasionally there are no individuals in a given group using the cut points. 
  group.exists<-paste0("strata(risk0_group)=risk0_group=",k)%in%summary(km.0.grp,cens=T)$strata
  
  if(group.exists){
    step.grp=stepfun(km.0.grp$time[summary(km.0.grp,cens=T)$strata==paste0("strata(risk0_group)=risk0_group=",k)],
                     c(1,km.0.grp$surv[summary(km.0.grp,cens=T)$strata==paste0("strata(risk0_group)=risk0_group=",k)]))
    risk0_true_grp[k]<-1-step.grp(5)
  }
}
calib_risk0_group[i,,3]<-risk0_true_grp

#-----
#true risk in each decile of risk - "always treated"

km.1.grp=survfit(Surv(T.A1,D.A1)~strata(risk1_group),data=dat.cf)

risk1_true_grp<-rep(NA,10)
for(k in 1:10){
  #This line is because occasionally there are no individuals in a given group using the cut points. 
  group.exists<-paste0("strata(risk1_group)=risk1_group=",k)%in%summary(km.1.grp,cens=T)$strata
  
  if(group.exists){
    step.grp=stepfun(km.1.grp$time[summary(km.1.grp,cens=T)$strata==paste0("strata(risk1_group)=risk1_group=",k)],
                     c(1,km.1.grp$surv[summary(km.1.grp,cens=T)$strata==paste0("strata(risk1_group)=risk1_group=",k)]))
    risk1_true_grp[k]<-1-step.grp(5)
  }
}
calib_risk1_group[i,,3]<-risk1_true_grp

#------------------------------------
#------------------------------------
#RESTRICTED SUBSET APPROACH
#Mean risk at times 1-5 in the subsets of individuals who are untreated/treated up to that time
#------------------------------------
#------------------------------------

#mean risks at times 1-5 in the subset of people who are untreated up to that time, under the "never treated" strategy
calib_risk0[i,,4]=sapply(1:5,FUN=function(x){mean(risk0_exp[subset0.matrix[,x]==1,x])})

#mean risks at times 1-5 in the subset of people who are treated up to that time, under the "always treated" strategy
calib_risk1[i,,4]=sapply(1:5,FUN=function(x){mean(risk1_exp[subset1.matrix[,x]==1,x])})

#------------------------------------
#------------------------------------
#RESTRICTED SUBSET APPROACH
#TRUE risks at times 1:5 in the subsets of individuals who are untreated/treated up to that time
#------------------------------------
#------------------------------------

calib_risk0[i,,5]=risk0_obs_subset
calib_risk1[i,,5]=risk1_obs_subset

#------------------------------------
#------------------------------------
#RESTRICTED SUBSET APPROACH
#Mean predicted and observed risks within categories of predicted risks under the two treatment strategies: "always treated", "never treated"
#in the subsets of individuals who are untreated/treated up to their end of follow-up
#------------------------------------
#------------------------------------

#-----
#generate expected risk groups
#USING DATA SET-SPECIFIC CUT POINTS

n0.subset<-dim(dat.val[subset0.5==1,])[1]
n1.subset<-dim(dat.val[subset1.5==1,])[1]

cut_points_0_subset=c(0,quantile(risk0_exp[subset0.5==1,5],probs=seq(0.1,0.9,0.1)),1)
cut_points_1_subset=c(0,quantile(risk1_exp[subset1.5==1,5],probs=seq(0.1,0.9,0.1)),1)

risk0_group_subset=cut(risk0_exp[subset0.5==1,5],breaks=cut_points_0_subset,include.lowest = T,labels = F)
dat.risk0.grp.subset=data.frame(id=dat.val$id[subset0.5==1],risk0_group_subset)
dat.0.merge.subset=merge(dat.risk0.grp.subset,dat.val[subset0.5==1,],by="id")

risk1_group_subset=cut(risk1_exp[subset1.5==1,5],breaks=cut_points_1_subset,include.lowest = T,labels = F)
dat.risk1.grp.subset=data.frame(id=dat.val$id[subset1.5==1],risk1_group_subset)
dat.1.merge.subset=merge(dat.risk1.grp.subset,dat.val[subset1.5==1,],by="id")

#-----
#mean predicted risk within each decile

calib_risk0_group[i,,4]<-sapply(FUN=function(x){mean(risk0_exp[subset0.5==1,5][risk0_group_subset==x])},1:10)
calib_risk1_group[i,,4]<-sapply(FUN=function(x){mean(risk1_exp[subset1.5==1,5][risk1_group_subset==x])},1:10)

#-----
#observed risk in each decile of risk - "never treated"

km.0.grp.subset=survfit(Surv(T.obs,D.obs)~strata(risk0_group_subset),data=dat.0.merge.subset)

risk0_obs_grp_subset<-rep(NA,10)
for(k in 1:10){
  #This line is because occasionally there are no individuals in a given group using the cut points. 
  group.exists<-paste0("strata(risk0_group_subset)=risk0_group_subset=",k)%in%summary(km.0.grp.subset,cens=T)$strata
  
  if(group.exists){
    step.grp=stepfun(km.0.grp.subset$time[summary(km.0.grp.subset,cens=T)$strata==paste0("strata(risk0_group_subset)=risk0_group_subset=",k)],
                     c(1,km.0.grp.subset$surv[summary(km.0.grp.subset,cens=T)$strata==paste0("strata(risk0_group_subset)=risk0_group_subset=",k)]))
    risk0_obs_grp_subset[k]<-1-step.grp(5)
  }
}
calib_risk0_group[i,,5]<-risk0_obs_grp_subset

#-----
#observed risk in each decile of risk - "always treated"

km.1.grp.subset=survfit(Surv(T.obs,D.obs)~strata(risk1_group_subset),data=dat.1.merge.subset)

risk1_obs_grp_subset<-rep(NA,10)
for(k in 1:10){
  #This line is because occasionally there are no individuals in a given group using the cut points. 
  group.exists<-paste0("strata(risk1_group_subset)=risk1_group_subset=",k)%in%summary(km.1.grp.subset,cens=T)$strata
  
  if(group.exists){
    step.grp=stepfun(km.1.grp.subset$time[summary(km.1.grp.subset,cens=T)$strata==paste0("strata(risk1_group_subset)=risk1_group_subset=",k)],
                     c(1,km.1.grp.subset$surv[summary(km.1.grp.subset,cens=T)$strata==paste0("strata(risk1_group_subset)=risk1_group_subset=",k)]))
    risk1_obs_grp_subset[k]<-1-step.grp(5)
  }
}
calib_risk1_group[i,,5]<-risk1_obs_grp_subset
