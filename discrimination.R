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

#-----------------------
# create model for censoring on all observed data, use to calculate weights for standard censoring

km.cens.stand <- survfit(Surv(T.obs,D.obs==0)~1,data=dat.val)
km.cens.step <- stepfun(km.cens.stand$time,c(1,km.cens.stand$surv))

#-----------------
#-----------------
# 'never treated' strategy
#-----------------
#-----------------

#-----------------------
# artificially censor the event time and event status in dat.val

dat.val$T.cens.0 <- ifelse(dat.val$A.0==1,0,
                        ifelse(dat.val$A.1==1,1,
                               ifelse(dat.val$A.2==1, 2,
                                      ifelse(dat.val$A.3==1, 3,
                                             ifelse(dat.val$A.4==1,4,dat.val$T.obs)))))
dat.val$D.cens.0 <- ifelse(dat.val$A.0==1,0,
                        ifelse(dat.val$A.1==1,0,
                               ifelse(dat.val$A.2==1, 0,
                                      ifelse(dat.val$A.3==1, 0,
                                             ifelse(dat.val$A.4==1,0,dat.val$D.obs)))))

#-----------------------
# time-dependent weights for artificial censoring (ipw0 and ipw1) have been made in 'analysis_x_model_validation.R' at visit time points (and at T.obs)
# for cindex and auct we need for each id weights at each event time that occurred in subjects who are not artificially censored

# unique event time points for 'never treated'
event.times.0 <- sort(unique(dat.long.val$time.stop[dat.long.val$in.dat.0==1 & dat.long.val$event==1]))
# we add the evaluation time point t=5 (we actually use 4.999) as we need it for calculation of CD AUCt and Brier score (at t=5)
event.times.0<-c(event.times.0,4.9999)

# select the rows from dat.long.val relevant for never treated scenario (using indicator in.dat.0) 
# split on event times relevant for the 'never treated' scenario (note that data are already split on visit times)
# note this step takes few seconds computation time (dependent on n)
dat.0.split <- survSplit( Surv(time,time.stop,event) ~., data = dat.long.val %>% filter(in.dat.0==1), cut = event.times.0)

# calculate weights for standard censoring at the end time points
dat.0.split$ipw.othercens <- 1/km.cens.step(dat.0.split$time.stop)
# combine these weights with the weights for artificial censoring
dat.0.split$ipw.comb <- dat.0.split$ipw0 * dat.0.split$ipw.othercens

#--------------
# construct weights matrix needed for cindex / auct 'never treated' 
# select the weights at event time points + make sure subjects who are censored before the first event time point are kept in the dataset
# put weights in wide format
dat.0.wide<- dat.0.split %>% 
  filter(time.stop %in% event.times.0 | (time.stop < min(event.times.0))) %>% 
  select(c("id","time.stop","ipw.comb")) %>% 
  spread(time.stop, ipw.comb)  

# subjects with censoring time before first event time should not add a row but not a column to the weightsmatrix (only columns for event time points are needed)
n.cens.before.first.event.0 <- sum(dat.0.split$time.stop < min(event.times.0))
dat.0.wide <- as.matrix(dat.0.wide[,-1:-(1+n.cens.before.first.event.0)]) #one additional column is deleted ("id")

# the above weightsmatrix only contains with at least A0=0. Expand it so that it has rows of NA for people in who are directly censored at t=0
wt_matrix0 <- matrix(nrow=n,ncol=ncol(dat.0.wide))
ids.0<-unique(dat.long.val$id[dat.long.val$A==0])
wt_matrix0[ids.0,] <- dat.0.wide

# remove last column from weightmatrix (weights at 4.9999) as this is not needed in the cindex calculation
wt_matrix0_eventsonly <- wt_matrix0[,-ncol(wt_matrix0)] 

#------------
# c-index for the 'never treated' strategy
#------------

#'real' c-index (using counterfactual data for the never treated strategy)
disc_cindex[i,1] <- c_index_ties(time=dat.cf$T.A0, status=dat.cf$D.A0, risk=risk0_exp[,5], tau=5)
#rcorr.cens(1-risk0_exp[,5],with(dat.cf,Surv(T.A0,D.A0)))[["C Index"]] #same

#unweighted c-index from artificially censored validation data for the never treated strategy (similar to Harrell's C) (not in the table/figures)
disc_cindex[i,2] <- c_index_ties(time=dat.val$T.cens.0, status=dat.val$D.cens.0, risk=risk0_exp[,5], tau=5)

#KM weighted c-index from artificially censored validation data for the never treated strategy (similar to Uno's c-index) (not in the table/figures) 
disc_cindex[i,3] <- concordance(Surv(dat.val$T.cens.0, dat.val$D.cens.0) ~ risk0_exp[,5],
                                dat.val,
                                reverse = TRUE,
                                timewt = "n/G2")$concordance

# time-dependent weighted c-index from artificially censored validation data for the never treated strategy (our proposal)
disc_cindex[i,4] <- c_index_ties(time=dat.val$T.cens.0, status=dat.val$D.cens.0, risk=risk0_exp[,5], tau=5, weightmatrix = wt_matrix0_eventsonly)

#subset method using weights (estimated in the subset) to account for standard censoring
#the subgroup indicator "subset0.5" was made in the file in analysis_x_model_validation.R
disc_cindex[i,9] <- concordance(Surv(T.obs, D.obs) ~ risk0_exp[subset0.5==1,5],
            dat.val[subset0.5==1,],
            reverse = TRUE,
            timewt = "n/G2")$concordance

#------------
# AUCt for the 'never treated' strategy
#------------

#'real' C/D AUCt (using counterfactual data for the never treated strategy)
rCDauct0 <- wCD_AUCt(time=dat.cf$T.A0, status=dat.cf$D.A0, risk=risk0_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,1] <- rCDauct0$AUCt$AUC[rCDauct0$AUCt$time==4.9999]

#unweighted C/D AUCt from artificially censored validation data for the never treated strategy (not in the table/figures)
uwCDauct0 <- wCD_AUCt(time=dat.val$T.cens.0, status=dat.val$D.cens.0, risk=risk0_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,2] <- uwCDauct0$AUCt$AUC[uwCDauct0$AUCt$time==4.9999]

#KM weighted C/D AUCt from artificially censored validation data for the never treated strategy (similar to Uno's AUC)(not in the table/figures)
KMwCDauct0 <- timeROC(
  T = dat.val$T.cens.0,
  delta = dat.val$D.cens.0,
  marker = risk0_exp[,5],
  cause = 1,
  weighting = "marginal",
  times = 4.9999,
  iid = FALSE
)
disc_auct[i,3] <- KMwCDauct0$AUC[2]

#weighted C/D AUCt from artificially censored validation data for the never treated strategy (our proposal)
wCDauct0 <- wCD_AUCt(time=dat.val$T.cens.0, status=dat.val$D.cens.0, risk=risk0_exp[,5], plot = FALSE, weightmatrix = wt_matrix0, seq.time=4.9999)
disc_auct[i,4] <- wCDauct0$AUCt$AUC[wCDauct0$AUCt$time==4.9999]

#C/D AUCt based on subset who did not start treatment at any time point
#subset when using weights (estimated in subset) to account for standard censoring
KMwCDauct0.subset <- timeROC(
  T = dat.val$T.obs[subset0.5==1],
  delta = dat.val$D.obs[subset0.5==1],
  marker = risk0_exp[subset0.5==1,5],
  cause = 1,
  weighting = "marginal",
  times = 4.9999,
  iid = FALSE
)
disc_auct[i,9] <- KMwCDauct0.subset$AUC[2]

#------------
# Brier and scaled Brier (IPA) for the 'never treated' strategy
#------------

#'real' Brier and scaled Brier (IPA) using counterfactual data for the never treated strategy 
brier_raw[i,1] <- Brier(dat.cf$T.A0, dat.cf$D.A0, risk0_exp[,5], seq.time=4.9999)
brier_ipa[i,1] <- ipa(dat.cf$T.A0, dat.cf$D.A0, risk0_exp[,5], seq.time=4.9999)

# KM-weighted Brier and scaled Brier (IPA)  from artificial censored validation data for the never treated strategy (not in table/figures)
# note these weights only account for artificial censoring, not for standard censoring
sf <- survfit(Surv(T.cens.0,1-D.cens.0)~1, data=dat.val)
probvals <- unique(sf$surv)[1:5] #probs of remaining uncensored in [0,1), [1,2), [2,3), [3-4), [4,5)
dat.val$prob_still_in_fu <- ifelse(dat.val$T.cens.0<1, probvals[1],
                               ifelse(dat.val$T.cens.0<2, probvals[2],
                                      ifelse(dat.val$T.cens.0<3, probvals[3],
                                             ifelse(dat.val$T.cens.0<4, probvals[4], probvals[5]))))
brier_raw[i,2] <- Brier(dat.val$T.cens.0, dat.val$D.cens.0, risk0_exp[,5], seq.time=4.9999, weights=1/dat.val$prob_still_in_fu)
brier_ipa[i,2] <- ipa(dat.val$T.cens.0, dat.val$D.cens.0, risk0_exp[,5], seq.time=4.9999, weights=1/dat.val$prob_still_in_fu)

#td-weighted Brier and scaled Brier (IPA) from artificial censored data for the never treated strategy (our proposal)
# first derive a vector of weights 
# for individuals with event after art cens (T.cens.0==1) - weight at their own event time
# for individuals known without event by tau after art cens (T.cens.0==5) - weight at tau
# weights for others are not needed (also does not mind if NA or other number)
# taking the last ipw0 weight recorded in dat.long.val for each id suffices here as max time is tau in the data
# calculate last ipw0 value in dat.long.val
dat.long.val <- dat.long.val %>% group_by(id) %>% mutate(ipw0.last=last(ipw0)) 
#  add this value to dat.val
dat.val <- left_join(dat.val, dat.long.val %>% 
                         select(c("id","ipw0.last")) %>% 
                         filter(row_number()==1), 
                       by="id")

# combine with standard censoring
dat.val$ipw.othercens <- 1/km.cens.step(dat.val$T.obs)
dat.val$ipw.othercens[dat.val$T.obs==5] <- 1/km.cens.step(4.9999)
dat.val$ipw0.comb <- dat.val$ipw0.last * dat.val$ipw.othercens

brier_raw[i,3] <- Brier(dat.val$T.cens.0, dat.val$D.cens.0, risk=risk0_exp[,5], seq.time=4.9999, weights=dat.val$ipw0.comb)
brier_ipa[i,3] <- ipa(dat.val$T.cens.0, dat.val$D.cens.0, risk=risk0_exp[,5], seq.time=4.9999, weights=dat.val$ipw0.comb)

#Brier and scaled Brier (IPA) based on subset who did not receive treatment at any time point
#use weights to account for standard uninformative censoring in the subset
km.cens.stand.subset0 <- survfit(Surv(T.obs,D.obs==0)~1,data=dat.val[subset0.5==1,])
pred.uncensored.subset0 <- evalstep(km.cens.stand.subset0$time,km.cens.stand.subset0$surv,dat.val$T.obs[subset0.5==1])
pred.uncensored.subset0[dat.val$T.obs[subset0.5==1]==5] <- evalstep(km.cens.stand.subset0$time,km.cens.stand.subset0$surv,4.9999)
wt.stand.subset0 <- 1/pred.uncensored.subset0

brier_raw[i,7] <- Brier(dat.val$T.obs[subset0.5==1], dat.val$D.obs[subset0.5==1], risk0_exp[subset0.5==1,5], seq.time=4.9999, weights=wt.stand.subset0)
brier_ipa[i,7] <- ipa(dat.val$T.obs[subset0.5==1], dat.val$D.obs[subset0.5==1], risk0_exp[subset0.5==1,5], seq.time=4.9999, weights=wt.stand.subset0)

#-----------------
#-----------------
# risk under 'always treated' strategy
#-----------------
#-----------------

#-----------------
#artificially censor the time and event status in dat.val

dat.val$T.cens.1 <- ifelse(dat.val$A.0==0,0,
                       ifelse(dat.val$A.1==0,pmin(dat.val$T.obs,1),
                              ifelse(dat.val$A.2==0, pmin(dat.val$T.obs,2),
                                     ifelse(dat.val$A.3==0, pmin(dat.val$T.obs,3),
                                            ifelse(dat.val$A.4==0,pmin(dat.val$T.obs,4),dat.val$T.obs)))))
dat.val$D.cens.1 <- ifelse(dat.val$A.0==0,0,
                       ifelse(dat.val$A.1==0 & dat.val$T.cens.1==1,0,
                              ifelse(dat.val$A.2==0 & dat.val$T.cens.1==2, 0,
                                     ifelse(dat.val$A.3==0 & dat.val$T.cens.1==3, 0,
                                            ifelse(dat.val$A.4==0 & dat.val$T.cens.1==4,0,dat.val$D.obs)))))

#-----------------------
# weights for artificial censoring (ipw1) have been made in 'analysis_x_model_validation.R' at visit time points (and at T.obs)
# for cindex and auct we need for each id weights at each event time that occurred in subjects who are not artificially censored

# unique event time points for 'always treated'
event.times.1 <- sort(unique(dat.long.val$time.stop[dat.long.val$in.dat.1==1 & dat.long.val$event==1]))
# we add the evaluation time point t=5 (we actually use 4.999) as we need it for calculation of CD AUCt and Brier score (at t=5)
event.times.1<-c(event.times.1,4.9999)

# select the rows from dat.long.val relevant for always treated scenario (using indicator in.dat.1) 
# split on event times relevant for the 'always treated' scenario (note that data are already split on visit times)
# note this step takes few seconds computation time (dependent on n)
dat.1.split <- survSplit( Surv(time,time.stop,event) ~., data = dat.long.val %>% filter(in.dat.1==1), cut = event.times.1)

# calculate weights for standard censoring at the end time points
dat.1.split$ipw.othercens <- 1/km.cens.step(dat.1.split$time.stop)
# combine these weights with the weights for artificial censoring
dat.1.split$ipw.comb <- dat.1.split$ipw1 * dat.1.split$ipw.othercens

#--------------
# construct weights matrix needed for cindex / auct 'never treated' 
# select the weights at event time points + make sure subjects who are censored before the first event time point are kept in the dataset
# put weights in wide format
dat.1.wide<- dat.1.split %>% 
  filter(time.stop %in% event.times.1 | (time.stop < min(event.times.1))) %>% 
  select(c("id","time.stop","ipw.comb")) %>% 
  spread(time.stop, ipw.comb)  

# subjects with censoring time before first event time should not add a column to the weightsmatrix (only weights at event time points are needed)
n.cens.before.first.event.1 <- sum(dat.1.split$time.stop<min(event.times.1))
dat.1.wide <- as.matrix(dat.1.wide[,-1:-(1+n.cens.before.first.event.1)]) #one additional column is deleted ("id")

# the above weights matrix is only for people with at least A0=0. Expand it so that it has rows of NA for people in who are directly censored at t=0
wt_matrix1 <- matrix(nrow=n,ncol=ncol(dat.1.wide))
ids.1 <- dat.long.val %>% filter(A==1 & time==0) %>% pull("id")
wt_matrix1[ids.1,] <- dat.1.wide

# remove last column from weightmatrix (weights at 4.9999) as this is not needed in the cindex calculation
wt_matrix1_eventsonly <- wt_matrix1[,-ncol(wt_matrix1)] 

#------------
# c-index for the 'always treated' strategy
#------------

#'real' c-index using counterfactual data for the always treated strategy
disc_cindex[i,5] <- c_index_ties(time=dat.cf$T.A1, status=dat.cf$D.A1, risk=risk1_exp[,5], tau=5)
#rcorr.cens(1-risk1_exp[,5],with(dat,Surv(T.A1,D.A1)))[["C Index"]] #same

#unweighted c-index from artificially censored validation data for the always treated strategy (Harrell's C) (not in the table/figures)
disc_cindex[i,6] <- c_index_ties(time=dat.val$T.cens.1, status=dat.val$D.cens.1, risk=risk1_exp[,5], tau=5)

# KM weighted c-index from artificially censored validation data for the always treated strategy (Uno's C) (not in the table/figures)
disc_cindex[i,7] <- concordance(Surv(dat.val$T.cens.1, dat.val$D.cens.1) ~ risk1_exp[,5],
                                dat.val,
                                reverse = TRUE,
                                timewt = "n/G2")$concordance

#weighted c-index from artificially censored validation data for the always treated strategy (our proposal)  
disc_cindex[i,8] <- c_index_ties(time=dat.val$T.cens.1, status=dat.val$D.cens.1, risk=risk1_exp[,5], tau=5, weightmatrix=wt_matrix1_eventsonly)

#cindex for subset using weights (estimated in the subset) to account for standard censoring
#the indicator subset1.5 was created in analysis_x_model_validation.R
disc_cindex[i,10] <-concordance(Surv(T.obs, D.obs) ~ risk1_exp[subset1.5==1,5],
                                dat.val[subset1.5==1,],
                                reverse = TRUE,
                                timewt = "n/G2")$concordance

#------------
# AUCt for the 'always treated' strategy
#------------

#'real' C/D AUCt using counterfactual data for the always treated strategy
rCDauct1 <- wCD_AUCt(time=dat.cf$T.A1, status=dat.cf$D.A1, risk=risk1_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,5] <- rCDauct1$AUCt$AUC[rCDauct1$AUCt$time==4.9999]

#unweighted C/D AUCt from artificially censored validation data for the always treated strategy (not in the table/figures)
uwCDauct1 <- wCD_AUCt(time=dat.val$T.cens.1, status=dat.val$D.cens.1, risk=risk1_exp[,5], plot = FALSE, seq.time=4.9999)
disc_auct[i,6] <- uwCDauct1$AUCt$AUC[uwCDauct1$AUCt$time==4.9999]

#KM weighted C/D AUCt from artificially censored validation data for the always treated strategy (similar to Uno's AUC) (not in the table/figures)
KMwCDauct1 <- timeROC(
  T = dat.val$T.cens.1,
  delta = dat.val$D.cens.1,
  marker = risk1_exp[,5],
  cause = 1,
  weighting = "marginal",
  #times = seq(0.1,4.99,.1),
  times = 4.9999,
  iid = FALSE
)
disc_auct[i,7] <- KMwCDauct1$AUC[2]

#weighted C/D AUCt from artificially censored validation data for the always treated strategy (our proposal) 
wCDauct1 <- wCD_AUCt(time=dat.val$T.cens.1, status=dat.val$D.cens.1, risk=risk1_exp[,5], plot = FALSE, weightmatrix = wt_matrix1, seq.time=4.9999)
disc_auct[i,8] <- wCDauct1$AUCt$AUC[wCDauct1$AUCt$time==4.9999]

#C/D AUCt in subset when using weights (estimated in subset) to account for standard censoring
KMwCDauct1.subset <- timeROC(
  T = dat.val$T.obs[subset1.5==1],
  delta = dat.val$D.obs[subset1.5==1],
  marker = risk1_exp[subset1.5==1,5],
  cause = 1,
  weighting = "marginal",
  times = 4.9999,
  iid = FALSE
)
disc_auct[i,10] <- KMwCDauct1.subset$AUC[2]

#------------
# Brier and scaled Brier (IPA) for the 'always treated' strategy
#------------

#'real' Brier and scaled Brier (IPA) from counterfactual data for the always treated strategy 
brier_raw[i,4] <- Brier(dat.cf$T.A1, dat.cf$D.A1, risk1_exp[,5], seq.time=4.9999)
brier_ipa[i,4] <- ipa(dat.cf$T.A1, dat.cf$D.A1, risk1_exp[,5], seq.time=4.9999)

#KM-weighted Brier and scaled Brier (IPA) from artificial censored data for the always treated strategy (not in table/figures) 
sf1 <- survfit(Surv(T.cens.1,1-D.cens.1)~1, data=dat.val)
#(only weighted for censoring at t=0, all uncensored observations get same weight)
weights.cens.1.KM <- 1/sf1$surv[1] #probs of remaining untreated at t=0
brier_raw[i,5] <- Brier(dat.val$T.cens.1, dat.val$D.cens.1, risk1_exp[,5], seq.time=4.9999, weights=rep(weights.cens.1.KM,n))
brier_ipa[i,5] <- ipa(dat.val$T.cens.1, dat.val$D.cens.1, risk1_exp[,5], seq.time=4.9999, weights=rep(weights.cens.1.KM,n))

#td-weighted Brier and scaled Brier (IPA) from artificial censored data for the never treated strategy (our proposal)
# first derive a vector of weights 
# for individuals with event after art cens (T.cens.0==1) - weight at their own event time
# for individuals known without event by tau after art cens (T.cens.0==5) - weight at tau
# weights for others are not needed (also does not mind if NA or other number)
# taking the last ipw1 weight recorded in dat.long.val for each id suffices here as max time is tau in the data
# calculate last ipw1 value in dat.long.val
dat.long.val <- dat.long.val %>% group_by(id) %>% mutate(ipw1.last=last(ipw1)) 
#  add this value to dat.val
dat.val <- left_join(dat.val, dat.long.val %>% 
                       select(c("id","ipw1.last")) %>% 
                       filter(row_number()==1), 
                     by="id")

# combine with standard censoring
dat.val$ipw1.comb <- dat.val$ipw1.last * dat.val$ipw.othercens

brier_raw[i,6] <- Brier(dat.val$T.cens.1, dat.val$D.cens.1, risk1_exp[,5], seq.time=4.9999, weights=dat.val$ipw1.comb)
brier_ipa[i,6] <- ipa(dat.val$T.cens.1, dat.val$D.cens.1, risk1_exp[,5], seq.time=4.9999, weights=dat.val$ipw1.comb)

#Brier and scaled Brier (IPA) based on subset who did receive treatment from time zero on
#use weights to account for standard uninformative censoring in the subset
km.cens.stand.subset1 <- survfit(Surv(T.obs,D.obs==0)~1,data=dat.val[subset1.5==1,])
pred.uncensored.subset1 <- evalstep(km.cens.stand.subset1$time,km.cens.stand.subset1$surv,dat.val$T.obs[subset1.5==1])
pred.uncensored.subset1[dat.val$T.obs[subset1.5==1]==5] <- evalstep(km.cens.stand.subset1$time,km.cens.stand.subset1$surv,4.9999)
wt.stand.subset1 <- 1/pred.uncensored.subset1

brier_raw[i,8] <- Brier(dat.val$T.obs[subset1.5==1], dat.val$D.obs[subset1.5==1], risk1_exp[subset1.5==1,5], seq.time=4.9999, weights=wt.stand.subset1)
brier_ipa[i,8] <- ipa(dat.val$T.obs[subset1.5==1], dat.val$D.obs[subset1.5==1], risk1_exp[subset1.5==1,5], seq.time=4.9999, weights=wt.stand.subset1)
