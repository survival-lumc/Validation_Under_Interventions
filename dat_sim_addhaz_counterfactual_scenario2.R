#------------------------------
#------------------------------
#Simulates counterfactual longitudinal data on covariates L and event times under 'always treated' and 'never treated, using an additive hazards model. 
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
alpha.0=0.3
alpha.A=-0.04
alpha.L=0.01
alpha.U=0.01
tfac=0.2

#---------------------
#---------------------
#generate event times and event indicators had each person been NEVER TREATED
#---------------------
#---------------------

A=matrix(nrow=n,ncol=n.visit)
L=matrix(nrow=n,ncol=n.visit)

L[,1]=dat.val$L.0
A[,1]=0
for(k in 2:n.visit){
  L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,sd.L)
  A[,k]=0
}

neghaz=0
T.A0=rep(NA,n)
for(k in 1:n.visit){
  u.t=runif(n,0,1)
  haz=alpha.0+alpha.A*A[,k]+alpha.L*(1-(k-1)*tfac)*L[,k]+alpha.U*U
  new.t=-log(u.t)/haz
  T.A0=ifelse(is.na(T.A0) & new.t<1 & haz>0,k-1+new.t,T.A0)
  neghaz=ifelse(sum(haz<0)>0,neghaz+sum(haz<0),0)
}
D.A0=ifelse(is.na(T.A0),0,1)
T.A0=ifelse(is.na(T.A0),5,T.A0)
neghaz

#---------------------
#---------------------
#generate event times and event indicators had each person been ALWAYS TREATED
#---------------------
#---------------------

A=matrix(nrow=n,ncol=n.visit)
L=matrix(nrow=n,ncol=n.visit)

L[,1]=dat.val$L.0
A[,1]=1
for(k in 2:n.visit){
  L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,sd.L)
  A[,k]=1
}

neghaz=0
T.A1=rep(NA,n)
for(k in 1:n.visit){
  u.t=runif(n,0,1)
  haz=alpha.0+alpha.A*A[,k]+alpha.L*(1-(k-1)*tfac)*L[,k]+alpha.U*U
  new.t=-log(u.t)/haz
  T.A1=ifelse(is.na(T.A1) & new.t<1 & haz>0,k-1+new.t,T.A1)
  neghaz=ifelse(sum(haz<0)>0,neghaz+sum(haz<0),0)
}
D.A1=ifelse(is.na(T.A1),0,1)
T.A1=ifelse(is.na(T.A1),5,T.A1)
neghaz

#-----
#Create data frame

dat.cf=data.frame(id=1:n,T.A0,D.A0,T.A1,D.A1,L.0=dat.val$L.0)