#---------------
#function for c-index with inverse prob weight option, allowing tied event times

c_index_ties <- function(time, status, risk, weightmatrix=NULL, tau)
  #Input:
  # time: event/censoring times
  # status: event indicators
  # risk: vector of risk probabilities by time horizon (typically tau) for each individual.
  # weightmatrix: matrix with ipc weights. Rows are subjects, columns are unique event time points (ordered)
  # tau: truncation time point, calculate c-index from zero up to (and including) tau
  #Output: 
  # c-index
{
  tt <- sort(unique(time[status == 1 & time <= tau]))
  if(is.null(weightmatrix))
  {
    weightmatrix=matrix(1,nrow=length(time),ncol = length(tt))
  }
  nt <- length(tt)            #number of unique event time points
  x <- risk 
  numsum <- denomsum <- 0
  for (i in 1:nt)                   #loop over unique event time points
  {
    ti <- tt[i]
    n1 <- intersect(which(time==ti) ,which(status==1)) #indices of cases at this time point
    n0 <- union(which(time>ti), intersect(which(time==ti),which(status==0)))  	#indices controls at this time point (patients with event time >ti plus patients censored at ti)
    nn1 <- length(n1)                 #number of cases
    nn0 <- length(n0)                 #number of controls
    x <- risk
    for (k in 1:nn1)                  #loop over the cases
    {
      xi <- x[n1][k] # risk of the k'th case at time point ti 
      numsum  <- numsum + weightmatrix[n1,i][k] * (x[n0]<xi)%*%weightmatrix[n0,i] + 0.5 * weightmatrix[n1,i][k] * (x[n0] == xi) %*% weightmatrix[n0,i]
      denomsum <- denomsum + weightmatrix[n1,i][k]*sum(weightmatrix[n0,i]) 
      
    }
    
  }		
  cindex <- numsum/denomsum
  return(cindex)
}


#---------------
#function for cumulative dynamic AUCt with ipc weights

wCD_AUCt <- function (time, status, risk, seq.time, plot = T, weightmatrix=NULL, xlim=5) 
  #Input:
  # time: event/censoring times
  # status: event indicators
  # risk: risk is a vector of risks for each subject (most implementations of C/D AUC just allow a vector, eg timeROC). But could be a matrix with predicitons by each unique event time point (our risk_obs_allt does not contain these for all observed event time points, only up to last follow up of each patient)
  # weightsmatrix: rows are subjects, weights are ordered unique event time points plus the evaluation time point appended (t=4.99 in this case, so weight at t=4)
  # seq.time: time vector where you want to calculate C/D AUC (must be ordered and smallest time should be equal or larger than first event time)
  # xlim: plotting parameter
  #Output: 
  # C/D AUCt values at seq.times
  # plot
{
  tt <- sort(unique(time[status==1])) #unique event time points
  if(is.null(weightmatrix))
  {
    weightmatrix=matrix(1,nrow=length(time),ncol = length(tt))
  }
  nseq <- length(seq.time)
  AUCt <- rep(NA,nseq)                  #vector to save AUCt in
  for (i in 1:nseq)                   #loop over time points where you want to calculate C/D AUCt
  {
    numsum <- denomsum <- 0
    ti <- seq.time[i]                 #ith unique time point where you want to calcuate C/D AUCt
    #tti <- which(tt<=ti))         #not used could be index among tt of this seq.time (then tt should be expanded by seq.time) and used in the weights of the controls below
    n1 <- intersect(which(time<=ti) ,which(status==1)) #indices of cumulative(!) cases at this time point
    n0 <- union(which(time>ti), intersect(which(time==ti),which(status==0)))  	#indices of controls at this time point (patients with event time >ti plus patients censored at ti)
    nn1 <- length(n1)                 #number of cases
    nn0 <- length(n0)                 #number of controls
    x <- risk                         #x is vector of risks to be evaluated (could be specific to i'th unique time point)
    for (k in 1:nn1)                  #loop over the cases
    {
      n1k <- n1[k]                    #index in original data of the k'th case
      xi <- x[n1k]                    # risk of the k'th case 
      ttk <- which(tt==time[n1k])     #index among unique event times for the k'th case
      numsum   <- numsum   + weightmatrix[n1k,ttk] * (x[n0]<xi)%*%weightmatrix[n0,ncol(weightmatrix)] + 0.5 * weightmatrix[n1k,ttk] * (x[n0] == xi) %*% weightmatrix[n0,ncol(weightmatrix)]
      denomsum <- denomsum + weightmatrix[n1k,ttk] * sum(weightmatrix[n0,ncol(weightmatrix)]) 
      #note that the weight for the case is evaluated at it's own unique event time point
      #the weight for the controls are evaluated at the time point where you want to calculate C/D AUCt
    }
    AUCt[i] <- numsum/denomsum
  }		
  
  if(plot) 
  {
    plot(seq.time,AUCt,xlab = "Time t", ylab = "C/D AUC(t)", type='l', xlim =c(0,xlim), ylim=c(0.5,1), main="cumulative dynamic AUCt")
    abline(h=0.5,lty=3)
  }
  return(list(AUCt=data.frame(time=seq.time,AUC=AUCt)))
}


