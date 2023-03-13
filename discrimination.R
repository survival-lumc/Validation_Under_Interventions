#----------------------------------
#----------------------------------
#In the current file we obtain estimates of c-index, AUCt, Brier score and weighted Brier score
#Output are four matrices containing
#
#disc_cindex
#   1. 'real' c-index using counterfactual data for the never treated strategy
#   2. unweighted c-index from artificially censored validation data for the never treated strategy (similar to Harrell's C) (not in the table/figures)
#   3. KM weighted c-index from artificially censored validation data for the never treated strategy (similar to Uno's c) (not in the table/figures)
#   4. time-dependent weighted c-index from artificially validation censored data for the never treated strategy (our proposal)
#   5. 'real' c-index using counterfactual data for the always treated strategy
#   6. unweighted c-index from artificially censored validation data for the always treated strategy (Harrell's C) (not in the table/figures)
#   7. KM weighted c-index from artificially censored validation data for the always treated strategy (Uno's C) (not in the table/figures)
#   8. time-dependent weighted c-index from artificially censored validation data for the always treated strategy (our proposal) 
#   9. c-index based on subset who did not start treatment at any time point
#   10.c-index based on subset who did receive treatment from time 0 onwards
#
#disc_auct
#   1  'real' C/D AUCt using counterfactual data for the never treated strategy
#   2. unweighted C/D AUCt from artificially censored validation data for the never treated strategy (not in the table/figures)
#   3. KM weighted C/D AUCt from artificially censored validation data for the never treated strategy (similar to Uno's AUC)(not in the table/figures)
#   4. time-dependent weighted C/D AUCt from artificially validation censored data for the never treated strategy (our proposal)
#   5. 'real' C/D AUCt using counterfactual data for the always treated strategy
#   6. unweighted C/D AUCt from artificially censored validation data for the always treated strategy (not in the table/figures)
#   7. KM weighted C/D AUCt from artificially censored validation data for the always treated strategy (similar to Uno's AUC) (not in the table/figures)
#   8. time-dependent weighted C/D AUCt from artificially censored validation data for the always treated strategy (our proposal) 
#   9. C/D AUCt based on subset who did not start treatment at any time point
#   10.C/D AUCt based on subset who did receive treatment from time 0 onwards
#
#brier_raw
#   1. 'real' Brier using counterfactual data for the never treated strategy 
#   2. KM-weighted Brier from artificial censored validation data for the never treated strategy (not in table/figures)
#   3. td-weighted Brier from artificial censored data for the never treated strategy (our proposal)
#   4. 'real' Brier from counterfactual data for the always treated strategy 
#   5. KM-weighted Brier from artificial censored data for the always treated strategy (not in table/figures)
#   6. weighted (for censoring at t=0 dependent on L0) Brier from artificial censored data for the always treated strategy (our proposal)
#   7. Brier based on subset who did not receive treatment at any time point
#   8. Brier based on subset who did receive treatment from time zero on
#
#brier_ipa same as brier_raw but now for Brier % improvement compared to null model
#-----------------
#-----------------

#-----------------
#-----------------
# risk under 'never treated' strategy
#-----------------
#-----------------

#-----------------------
# artificially censor the data

dat.val$T.obs.A0 <- ifelse(dat.val$A.0==1,0,
                       ifelse(dat.val$A.1==1,1,
                              ifelse(dat.val$A.2==1, 2,
                                     ifelse(dat.val$A.3==1, 3,
                                            ifelse(dat.val$A.4==1,4,dat.val$T.obs)))))
dat.val$D.obs.A0 <- ifelse(dat.val$A.0==1,0,
                       ifelse(dat.val$A.1==1,0,
                              ifelse(dat.val$A.2==1, 0,
                                     ifelse(dat.val$A.3==1, 0,
                                            ifelse(dat.val$A.4==1,0,dat.val$D.obs)))))

#-----------------------
# create time-dependent weights

#we need weights at each event time plus at the evaluation time point (t=5)
#given how data were simulated, weights are constant in between visit time, ie no treatment start in between),
#therefore we can use the weights at time point of last visit before the event time
#we add time point 4.9999 to the weightmatrix as we need it for calculation of CD AUCt (at t=5)

wt.mod = glm(A~L,family="binomial",data=dat.long.val[dat.long.val$Alag1==0,]) #weights model
pred.wt0 = predict(wt.mod,type = "response",newdata = dat.long.val) #prob of being artificially censored in certain time interval
dat.long.val$wt0 = (1-pred.wt0) #prob of not being artificially censored in certain time interval
dat.long.val$wt0 = ifelse(dat.long.val$A==1,NA,dat.long.val$wt0) #weights should not be used after tx start
dat.long.val$wt.cum0 = ave(dat.long.val$wt0,dat.long.val$id,FUN=cumprod) #prob of not being artificially censored cumulative

event.times.A0 <- with(dat.val, sort(T.obs.A0[D.obs.A0==1&T.obs.A0<=5])) #sorted event times after artificial censoring
event.times.A0 <- c(event.times.A0, 4.9999) #add time point 4.9999 for CD AUCt at time 4.9999

wt_dat <- data.frame(
  id = rep( sort(unique(dat.val$id)),     	    # replicate patient ID for each
            each = length(event.times.A0) ),    # unique event timepoints (plus 4.9999)
  event.time = rep((event.times.A0),
                   each = length(id)))

wt_dat$time=floor(wt_dat$event.time)          # add time point of last visit before the event time point to this matrix for merging with dat.long
dat.long.val.select <- select(dat.long.val,c("id","time","wt.cum0") )
wt_dat <- left_join(wt_dat, dat.long.val.select, by=c("id","time"))

wt_dat <- wt_dat[with(wt_dat, order(id, event.time)), ] #order according to event times again (order gets lost after merging)
wt_dat <- wt_dat[,-3] #remove time of visit from the data frame
wt_matrix <- pivot_wider(names_from = event.time, values_from = wt.cum0, data=wt_dat) #from long to wide
wt_matrix <- wt_matrix[,-1] #remove id from the matrix
wt_matrix0 <- 1/wt_matrix #weights are the inverse of prob of not being artificially censored
 
#------------
# c-index for the 'never treated' strategy
#------------

#'real' c-index (using counterfactual data for the never treated strategy)
disc_cindex[i,1] <- c_index_ties(time=dat.cf$T.A0, status=dat.cf$D.A0, risk=risk0_exp[,5], tau=5)
#rcorr.cens(1-risk0_exp[,5],with(dat.cf,Surv(T.A0,D.A0)))[["C Index"]] #same

#unweighted c-index from artificially censored validation data for the never treated strategy (similar to Harrell's C) (not in the table/figures)
disc_cindex[i,2] <- c_index_ties(time=dat.val$T.obs.A0, status=dat.val$D.obs.A0, risk=risk0_exp[,5], tau=5)

#KM weighted c-index from artificially censored validation data for the never treated strategy (similar to Uno's c-index) (not in the table/figures) 
disc_cindex[i,3] <- concordance(Surv(dat.val$T.obs.A0, dat.val$D.obs.A0) ~ risk0_exp[,5],
                                dat.val,
                                reverse = TRUE,
                                timewt = "n/G2")$concordance

# time-dependent weighted c-index from artificially validation censored data for the never treated strategy (our proposal)
wt_matrix0_eventsonly <- wt_matrix0[,-ncol(wt_matrix0)] #calculations for c-index only use weights at event times so last column of the weight-matrix is not needed
disc_cindex[i,4] <- c_index_ties(time=dat.val$T.obs.A0, status=dat.val$D.obs.A0, risk=risk0_exp[,5], tau=5, weightmatrix = wt_matrix0_eventsonly)

# subset method (using data for only those who never started treatment)
uncensored0 <- (dat.val$T.obs.A0>4.9999 | dat.val$D.obs.A0>0)
disc_cindex[i,9] <- c_index_ties(time=dat.val$T.obs[uncensored0], status=dat.val$D.obs[uncensored0], risk=risk0_exp[uncensored0,5], tau=5)
#rcorr.cens(1-risk0_exp[uncensored0,5],with(dat.val[uncensored0,],Surv(T.obs,D.obs)))[["C Index"]] #same

#------------
# AUCt for the 'never treated' strategy
#------------

#'real' C/D AUCt (using counterfactual data for the never treated strategy)
rCDauct0 <- wCD_AUCt(time=dat.cf$T.A0, status=dat.cf$D.A0, risk=risk0_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,1] <- rCDauct0$AUCt$AUC[rCDauct0$AUCt$time==4.9999]

#unweighted C/D AUCt from artificially censored validation data for the never treated strategy (similar to Harrell's C) (not in the table/figures)
uwCDauct0 <- wCD_AUCt(time=dat.val$T.obs.A0, status=dat.val$D.obs.A0, risk=risk0_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,2] <- uwCDauct0$AUCt$AUC[uwCDauct0$AUCt$time==4.9999]

#KM weighted C/D AUCt from artificially censored validation data for the never treated strategy (similar to Uno's AUC)(not in the table/figures)
KMwCDauct0 <- timeROC(
  T = dat.val$T.obs.A0,
  delta = dat.val$D.obs.A0,
  marker = risk0_exp[,5],
  cause = 1,
  weighting = "marginal",
  #times = seq(0.1,4.99,.1),
  times = 4.9999,
  iid = FALSE
)
disc_auct[i,3] <- KMwCDauct0$AUC[2]

#unweighted C/D AUCt from artificially censored validation data for the never treated strategy (not in the table/figures)
wCDauct0 <- wCD_AUCt(time=dat.val$T.obs.A0, status=dat.val$D.obs.A0, risk=risk0_exp[,5], plot = FALSE, weightmatrix = wt_matrix0, seq.time=4.9999)
disc_auct[i,4] <- wCDauct0$AUCt$AUC[wCDauct0$AUCt$time==4.9999]

#C/D AUCt based on subset who did not start treatment at any time point
restr_CDAUCt0 <- wCD_AUCt(time=dat.val$T.obs[uncensored0], status=dat.val$D.obs[uncensored0], risk=risk0_exp[uncensored0,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,9] <- restr_CDAUCt0$AUCt$AUC[restr_CDAUCt0$AUCt$time==4.9999]

#------------
# Brier and scaled Brier (IPA) for the 'never treated' strategy
#------------

#'real' Brier and scaled Brier (IPA) using counterfactual data for the never treated strategy 
brier_raw[i,1] <- Brier(dat.cf$T.A0, dat.cf$D.A0, risk0_exp[,5], seq.time=4.9999)
brier_ipa[i,1] <- ipa(dat.cf$T.A0, dat.cf$D.A0, risk0_exp[,5], seq.time=4.9999)

# KM-weighted Brier and scaled Brier (IPA)  from artificial censored validation data for the never treated strategy (not in table/figures)
sf <- survfit(Surv(T.obs.A0,1-D.obs.A0)~1, data=dat.val)
probvals <- unique(sf$surv)[1:5] #probs of remaining uncensored in [0,1), [1,2), [2,3), [3-4), [4,5)
dat.val$prob_still_in_fu <- ifelse(dat.val$T.obs.A0<1, probvals[1],
                               ifelse(dat.val$T.obs.A0<2, probvals[2],
                                      ifelse(dat.val$T.obs.A0<3, probvals[3],
                                             ifelse(dat.val$T.obs.A0<4, probvals[4], probvals[5]))))
brier_raw[i,2] <- Brier(dat.val$T.obs.A0, dat.val$D.obs.A0, risk0_exp[,5], seq.time=4.9999, weights=1/dat.val$prob_still_in_fu)
brier_ipa[i,2] <- ipa(dat.val$T.obs.A0, dat.val$D.obs.A0, risk0_exp[,5], seq.time=4.9999, weights=1/dat.val$prob_still_in_fu)

#td-weighted Brier and scaled Brier (IPA) from artificial censored data for the never treated strategy (our proposal)
weights.obs.A0 <- rep(NA,n)
tt <- sort(unique(dat.val$T.obs.A0[dat.val$D.obs.A0==1])) #unique event time points
n1 <- which(dat.val$D.obs.A0==1)              #indices of events
for (j in 1: length(tt))
{
  n1k <- n1[j]                    #index in original data of the k'th case
  ttk <- which(tt==dat.val$T.obs.A0[n1k])     #index in the vector of unique event times
  weights.obs.A0[n1k] <- wt_matrix0[n1k,ttk] # weight for n1kth patient, at their own event time
}
weights.obs.A0[dat.val$D.obs.A0==0] <- wt_matrix0[dat.val$D.obs.A0==0,ncol(wt_matrix0)] # censoring only occurs at t=5, take weight at that time point (censored observations before t=5 get weight=0)
brier_raw[i,3] <- Brier(dat.val$T.obs.A0, dat.val$D.obs.A0, risk=risk0_exp[,5], seq.time=4.9999, weights=weights.obs.A0)
brier_ipa[i,3] <- ipa(dat.val$T.obs.A0, dat.val$D.obs.A0, risk=risk0_exp[,5], seq.time=4.9999, weights=weights.obs.A0)

#Brier and scaled Brier (IPA) based on subset who did not receive treatment at any time point
brier_raw[i,7] <- Brier(dat.val$T.obs[uncensored0], dat.val$D.obs[uncensored0], risk0_exp[uncensored0,5], seq.time=4.9999)
brier_ipa[i,7] <- ipa(dat.val$T.obs[uncensored0], dat.val$D.obs[uncensored0], risk0_exp[uncensored0,5], seq.time=4.9999)

#-----------------
#-----------------
# risk under 'always treated' strategy
#-----------------
#-----------------

#-----------------
#artificially censor the data

dat.val$T.obs.A1 <- ifelse(dat.val$A.0==0,0,
                       ifelse(dat.val$A.1==0,pmin(dat.val$T.obs,1),
                              ifelse(dat.val$A.2==0, pmin(dat.val$T.obs,2),
                                     ifelse(dat.val$A.3==0, pmin(dat.val$T.obs,3),
                                            ifelse(dat.val$A.4==0,pmin(dat.val$T.obs,4),dat.val$T.obs)))))
dat.val$D.obs.A1 <- ifelse(dat.val$A.0==0,0,
                       ifelse(dat.val$A.1==0 & dat.val$T.obs.A1==1,0,
                              ifelse(dat.val$A.2==0 & dat.val$T.obs.A1==2, 0,
                                     ifelse(dat.val$A.3==0 & dat.val$T.obs.A1==3, 0,
                                            ifelse(dat.val$A.4==0 & dat.val$T.obs.A1==4,0,dat.val$D.obs)))))

#-----------------
#create weights for dependent censoring based on L0

wt.mod.1 <- glm(A~L,family="binomial",data=dat.long.val[dat.long.val$time==0,])
pred.wt1=predict(wt.mod.1,type = "response",newdata = dat.long.val[dat.long.val$time==0,]) #prob of not being artificially censored (getting treatment means not censored)
#wt1=(1-pred.wt1) #not needed here only at the untreated scenario 
weights.obs.A1 <- 1/pred.wt1 #probs of remaining uncensored at t=0
event.times.A1 <- with(dat.val, sort(T.obs.A1[D.obs.A1==1&T.obs.A1<=5]))
wt_matrix1 <- matrix(rep(weights.obs.A1, length(event.times.A1)),nrow=n,ncol = length(event.times.A1))

#------------
# c-index for the 'always treated' strategy
#------------

#'real' c-index using counterfactual data for the always treated strategy
disc_cindex[i,5] <- c_index_ties(time=dat.cf$T.A1, status=dat.cf$D.A1, risk=risk1_exp[,5], tau=5)
#rcorr.cens(1-risk1_exp[,5],with(dat,Surv(T.A1,D.A1)))[["C Index"]] #same

#unweighted c-index from artificially censored validation data for the always treated strategy (Harrell's C) (not in the table/figures)
disc_cindex[i,6] <- c_index_ties(time=dat.val$T.obs.A1, status=dat.val$D.obs.A1, risk=risk1_exp[,5], tau=5)

# KM weighted c-index from artificially censored validation data for the always treated strategy (Uno's C) (not in the table/figures)
disc_cindex[i,7] <- concordance(Surv(dat.val$T.obs.A1, dat.val$D.obs.A1) ~ risk1_exp[,5],
                                dat.val,
                                reverse = TRUE,
                                timewt = "n/G2")$concordance

#L0 weighted c-index from artificially censored validation data for the always treated strategy (our proposal)  
disc_cindex[i,8] <- c_index_ties(time=dat.val$T.obs.A1, status=dat.val$D.obs.A1, risk=risk1_exp[,5], tau=5, weightmatrix=wt_matrix1)

#c-index based on subset who did receive treatment from time 0 onwards
uncensored1 <- (dat.val$T.obs.A1>4.9999 | dat.val$D.obs.A1>0)
disc_cindex[i,10] <- c_index_ties(time=dat.val$T.obs[uncensored1], status=dat.val$D.obs[uncensored1], risk=risk1_exp[uncensored1,5], tau=5)
#rcorr.cens(1-risk1_exp[uncensored1,5],with(dat.val[uncensored1,],Surv(T.obs,D.obs)))[["C Index"]] #same

#------------
# AUCt for the 'always treated' strategy
#------------

#'real' C/D AUCt using counterfactual data for the always treated strategy
rCDauct1 <- wCD_AUCt(time=dat.cf$T.A1, status=dat.cf$D.A1, risk=risk1_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,5] <- rCDauct1$AUCt$AUC[rCDauct1$AUCt$time==4.9999]

#unweighted C/D AUCt from artificially censored validation data for the always treated strategy (not in the table/figures)
uwCDauct1 <- wCD_AUCt(time=dat.val$T.obs.A1, status=dat.val$D.obs.A1, risk=risk1_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,6] <- uwCDauct1$AUCt$AUC[uwCDauct1$AUCt$time==4.9999]

#KM weighted C/D AUCt from artificially censored validation data for the always treated strategy (similar to Uno's AUC) (not in the table/figures)
KMwCDauct1 <- timeROC(
  T = dat.val$T.obs.A1,
  delta = dat.val$D.obs.A1,
  marker = risk1_exp[,5],
  cause = 1,
  weighting = "marginal",
  #times = seq(0.1,4.99,.1),
  times = 4.9999,
  iid = FALSE
)
disc_auct[i,7] <- KMwCDauct1$AUC[2]

#L0 weighted C/D AUCt from artificially censored validation data for the always treated strategy (our proposal) 
wCDauct1 <- wCD_AUCt(time=dat.val$T.obs.A1, status=dat.val$D.obs.A1, risk=risk1_exp[,5], plot = FALSE, weightmatrix = wt_matrix1, seq.time=4.9999)
disc_auct[i,8] <- wCDauct1$AUCt$AUC[wCDauct1$AUCt$time==4.9999]

#C/D AUCt based on subset who did receive treatment from time 0 onwards
restr_CDAUCt1 <- wCD_AUCt(time=dat.val$T.obs[uncensored1], status=dat.val$D.obs[uncensored1], risk=risk1_exp[uncensored1,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,10] <- restr_CDAUCt1$AUCt$AUC[restr_CDAUCt1$AUCt$time==4.9999]

#------------
# Brier and scaled Brier (IPA) for the 'always treated' strategy
#------------

#'real' Brier and scaled Brier (IPA) from counterfactual data for the always treated strategy 
brier_raw[i,4] <- Brier(dat.cf$T.A1, dat.cf$D.A1, risk1_exp[,5], seq.time=4.9999)
brier_ipa[i,4] <- ipa(dat.cf$T.A1, dat.cf$D.A1, risk1_exp[,5], seq.time=4.9999)

#KM-weighted Brier and scaled Brier (IPA) from artificial censored data for the always treated strategy (not in table/figures) 
#only weighted for censoring at t=0, all uncensored observations get same weight)
sf1 <- survfit(Surv(T.obs.A1,1-D.obs.A1)~1, data=dat.val)
#'dropout' is only at t=0, after that all will remain untreated up to t=5
weights.obs.A1.KM <- 1/sf1$surv[1] #probs of remaining untreated at t=0
brier_raw[i,5] <- Brier(dat.val$T.obs.A1, dat.val$D.obs.A1, risk1_exp[,5], seq.time=4.9999, weights=rep(weights.obs.A1.KM,n))
brier_ipa[i,5] <- ipa(dat.val$T.obs.A1, dat.val$D.obs.A1, risk1_exp[,5], seq.time=4.9999, weights=rep(weights.obs.A1.KM,n))

#weighted (for censoring at t=0 dependent on L0) Brier and scaled Brier (IPA) from artificial censored data for the always treated strategy (our proposal)
brier_raw[i,6] <- Brier(dat.val$T.obs.A1, dat.val$D.obs.A1, risk1_exp[,5], seq.time=4.9999, weights=weights.obs.A1)
brier_ipa[i,6] <- ipa(dat.val$T.obs.A1, dat.val$D.obs.A1, risk1_exp[,5], seq.time=4.9999, weights=weights.obs.A1)

#Brier and scaled Brier (IPA) based on subset who did receive treatment from time zero on
brier_raw[i,8] <- Brier(dat.val$T.obs[uncensored1], dat.val$D.obs[uncensored1], risk1_exp[uncensored1,5], seq.time=4.9999)
brier_ipa[i,8] <- ipa(dat.val$T.obs[uncensored1], dat.val$D.obs[uncensored1], risk1_exp[uncensored1,5], seq.time=4.9999)
