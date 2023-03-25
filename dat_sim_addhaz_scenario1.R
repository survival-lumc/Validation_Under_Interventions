#------------------------------
#------------------------------
#simulates longitudinal data on treatment A and covariates L at 5 time points, and then generates event times according to an additive hazards model. 
#U is an individual frailty. 
#People who do not have the event are censored at time 5.
#
#SCENARIO 2: basis scenario
#------------------------------
#------------------------------

#----
#number of visits (K+1)

n.visit=5

#----
#expit function

expit=function(x){exp(x)/(1+exp(x))}

#------------------
#parameter values
#------------------
  
#SD(U)
sd.U<-2

#SD of L|U
sd.L=4

#Mean of L: U+mu.L
mu.L<-10

#model for A|L
gamma.0=-2
gamma.L=0.1

#model for hazard
alpha.0=0.2
alpha.A=-0.04
alpha.L=0.01
alpha.U=0.01
tfac=0.2

#------------------
#simulate data
#------------------

#----
#generate U, A, L

A=matrix(nrow=n,ncol=n.visit)
L=matrix(nrow=n,ncol=n.visit)

U=rnorm(n,0,sd.U)
L[,1]=rnorm(n,U+mu.L,sd.L)
A[,1]=rbinom(n,1,expit(gamma.0+gamma.L*L[,1]))
for(k in 2:n.visit){
  L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,sd.L)
  A[,k]=ifelse(A[,k-1]==1,1,rbinom(n,1,expit(gamma.0+gamma.L*L[,k])))
}

#----
#generate event times T.obs, and event indicators D.obs

T.obs=rep(NA,n)

neghaz=0
for(k in 1:n.visit){
  u.t=runif(n,0,1)
  haz=alpha.0+alpha.A*A[,k]+alpha.L*(1-(k-1)*tfac)*L[,k]+alpha.U*U
  new.t=-log(u.t)/haz
  T.obs=ifelse(is.na(T.obs) & new.t<1 & haz>0,k-1+new.t,T.obs)#the haz>0 is just used to deal with tiny possibility (under this data generating mechanism) the hazard could go negative. 
  neghaz=ifelse(sum(haz<0)>0,neghaz+sum(haz<0),0)
}
D.obs=ifelse(is.na(T.obs),0,1)
T.obs=ifelse(is.na(T.obs),5,T.obs)
neghaz

#----
#generate standard (uninformative) censoring times and adjust event times T.obs, and event indicator D.obs accordingly

if(censoring==TRUE)
{
  u.c=runif(n,0,1)
  haz.c=0.1 #without events about 40% censored before t=5
  #haz=0.05 #without events about 22% censored before t=5
  C =-log(u.c)/haz.c
  #sum(C<5)/n

  T.obs=ifelse(T.obs>C ,C, T.obs)
  D.obs=ifelse(T.obs<C & T.obs < 5,1,0)
}
  
#-----
#Create data frame

colnames(A)=paste0("A.",0:4)
colnames(L)=paste0("L.",0:4)
dat=data.frame(id=1:n,T.obs,D.obs,A,L)

#-----
#set A to 0 in time periods after event/censoring
 
dat$A.1=ifelse(dat$T.obs<1,0,dat$A.1)
dat$A.2=ifelse(dat$T.obs<2,0,dat$A.2)
dat$A.3=ifelse(dat$T.obs<3,0,dat$A.3)
dat$A.4=ifelse(dat$T.obs<4,0,dat$A.4)

#------------------
#some summaries: may be useful if you wish to change the parameter values used above, to consider other scenarios.
#------------------
 
#proportion always treated
always.treat=A[,1]+A[,2]+A[,3]+A[,4]+A[,5]

#proportion never treated
never.treat=(1-A[,1])+(1-A[,2])+(1-A[,3])+(1-A[,4])+(1-A[,5])

tabyl(always.treat)
tabyl(never.treat)
tabyl(dat$D.obs)

#------------------
#Reshape data into 'long' format (multiple rows per individual: 1 row for each visit)
#------------------

dat.long=reshape(data = dat,varying=c(paste0("A.",0:4),paste0("L.",0:4)),direction="long",idvar="id")
dat.long=dat.long[order(dat.long$id,dat.long$time),]

#generate start and stop times for each row
dat.long$time.stop=dat.long$time+1

dat.long=dat.long[dat.long$time<dat.long$T.obs,]

dat.long$time.stop=ifelse(dat.long$time.stop>dat.long$T.obs,dat.long$T.obs,dat.long$time.stop)

dat.long$event=ifelse(dat.long$time.stop==dat.long$T.obs & dat.long$D.obs==1,1,0)

#visit number
dat.long$visit=ave(rep(1,dim(dat.long)[1]),dat.long$id,FUN=cumsum)

#generate lagged A values
dat.long=dat.long %>%
  group_by(id) %>%
  mutate(Alag1 = lag(A,n=1),Alag2 = lag(A,n=2),Alag3 = lag(A,n=3),Alag4 = lag(A,n=4)) %>%
  mutate(Alag1=replace_na(Alag1,0),Alag2=replace_na(Alag2,0),Alag3=replace_na(Alag3,0),Alag4=replace_na(Alag4,0))

#generate lagged L values
dat.long=dat.long %>%
  group_by(id) %>%
  mutate(Llag1 = lag(L,n=1),Llag2 = lag(L,n=2),Llag3 = lag(L,n=3),Llag4 = lag(L,n=4)) %>%
  mutate(Llag1=replace_na(Llag1,0),Llag2=replace_na(Llag2,0),Llag3=replace_na(Llag3,0),Llag4=replace_na(Llag4,0))

#baseline L
dat.long=dat.long %>%
  group_by(id) %>%
  mutate(L.baseline = first(L))
