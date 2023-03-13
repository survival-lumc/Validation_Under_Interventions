#------------------
# save some descriptives for each simulated dataset
#------------------

#---
# percentage of patient*visit combinations where treatment was given (development data)
always.treat=dat.dev$A.0+dat.dev$A.1+dat.dev$A.2+dat.dev$A.3+dat.dev$A.4
descriptives_scenario[i,1] <- sum(always.treat)/(n.visit*n)
colnames(descriptives_scenario)[1] <- c("prop_treated_visits")

#---
# percentage of patients who were never treated (development data)
descriptives_scenario[i,2] <- tabyl(always.treat)[1,3]
colnames(descriptives_scenario)[2] <- c("prop_never_treated")

#---
# true risk distribution for never treated scenario in all and in uncensored data
source("analysis_cox_model_truerisks.R")

descriptives_scenario[i,3] <- mean(risk0_true[,5])
descriptives_scenario[i,4] <- sd(risk0_true[,5])

uncensored0 <- dat.val$T.obs.A0>4.9999 | dat.val$D.obs.A0>0
censored0 <- !uncensored0
descriptives_scenario[i,5] <- mean(risk0_true[uncensored0,5])
descriptives_scenario[i,6] <- sd(risk0_true[uncensored0,5])

colnames(descriptives_scenario)[3:6] <- c("mean_truerisk_nevertreat","sd_truerisk_nevertreat", "mean_truerisk_nevertreat_uncensored", "sd_truerisk_nevertreat_uncensored")

#---
# calculate theoretical discrimination index (asymptotic formula for binary events)

m1 <- descriptives_scenario[i,7] <- mean(risk0_true[dat.cf$D.A0==1,5])
s1<- descriptives_scenario[i,8] <- sd(risk0_true[dat.cf$D.A0==1,5])
m0 <- descriptives_scenario[i,9] <- mean(risk0_true[dat.cf$D.A0==0,5])
s0 <- descriptives_scenario[i,10] <- sd(risk0_true[dat.cf$D.A0==0,5])
descriptives_scenario[i,11] <- pnorm((m1-m0)/sqrt(s1^2+s0^2))

colnames(descriptives_scenario)[7:11] <- c("mean_truerisk_nevertreated_cases", "sd_truerisk_nevertreated_cases", "mean_truerisk_nevertreated_controls", "sd_truerisk_nevertreated_controls", "theor_c_stat")

#---
# calculate summaries on numbers of events in counterfactual data and validation data after artificial censoring

descriptives_scenario[i,12] <- sum(dat.val$D.obs.A0)
descriptives_scenario[i,13] <- sum(dat.val$D.obs.A1)
descriptives_scenario[i,14] <- sum(dat.cf$D.A0)
descriptives_scenario[i,15] <- sum(dat.cf$D.A1)

colnames(descriptives_scenario)[12:15] <- c("n_events_dat_val0", "n_events_dat_val1", "n_events_dat_cf0", "n_events_dat_cf1")

#-----------------------------------
# data for descriptive KM plots

sf.dev <- survfit(Surv(T.obs, D.obs)~1, data=dat.dev )
sf.cf1 <- survfit(Surv(T.A1, D.A1)~1, data=dat.cf)
sf.cf0 <- survfit(Surv(T.A0, D.A0)~1, data=dat.cf)
sf.val <- survfit(Surv(T.obs, D.obs)~1, data=dat.val)
s <- list("development" = sf.dev,
          "cf treated" = sf.cf1,
          "cf untreated" = sf.cf0,
          "validation"=sf.val)

grid <- seq(from = 0, to = 5, by = 0.1)
risk.dev <- evalstep(sf.dev$time,1-sf.dev$surv,grid)
risk.dev[1] <- 0
risk.cf1 <- evalstep(sf.cf1$time,1-sf.cf1$surv,grid)
risk.cf1[1] <- 0
risk.cf0 <- evalstep(sf.cf0$time,1-sf.cf0$surv,grid)
risk.cf0[1] <- 0
risk.val <- evalstep(sf.val$time,1-sf.val$surv,grid)
risk.val[1] <- 0

descr_KM_curves <- data.frame(
  sim = rep(i,4*length(grid)),
  time = rep(grid,4),
  type = rep(c("development", "validation","cf_treated", "cf_untreated"),each=length(grid)),
  risk = c(risk.dev, risk.val, risk.cf1, risk.cf0))
  
descr_KM_curves_long <- rbind(descr_KM_curves_long, descr_KM_curves)



