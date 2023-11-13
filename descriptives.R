#------------------
# save some descriptives for each simulated dataset
#------------------

#---
# percentage of patient*visit combinations where treatment was given (development data)
always.treat=dat.dev$A.0+dat.dev$A.1+dat.dev$A.2+dat.dev$A.3+dat.dev$A.4
descriptives_scenario[i,1] <- sum(always.treat)/(n.visit*n)
colnames(descriptives_scenario)[1] <- c("prop_treated_visits_dev")

#---
# percentage of patients who were always treated (development data)
descriptives_scenario[i,2] <- tabyl(always.treat)[6,3]
colnames(descriptives_scenario)[2] <- c("prop_always_treated_dev")

#---
# percentage of patients who were never treated (development data)
descriptives_scenario[i,3] <- tabyl(always.treat)[1,3]
colnames(descriptives_scenario)[3] <- c("prop_never_treated_dev")

#---
# calculate summaries on numbers of events in counterfactual data and validation data after artificial censoring

descriptives_scenario[i,13] <- sum(dat.val$D.cens.0)
descriptives_scenario[i,14] <- sum(dat.val$D.cens.1)
descriptives_scenario[i,15] <- sum(dat.cf$D.A0)
descriptives_scenario[i,16] <- sum(dat.cf$D.A1)

colnames(descriptives_scenario)[13:16] <- c("n_events_dat_val0", "n_events_dat_val1", "n_events_dat_cf0", "n_events_dat_cf1")

#--- 
# calculate summaries on number of patients still at risk in last time interval (T>4) in counterfactual data and validation data after censoring
# so patients untreated at all 5 visits and no event before T=4

descriptives_scenario[i,17] <- sum(dat.val$T.cens.0>4)
descriptives_scenario[i,18] <- sum(dat.val$T.cens.1>4)
descriptives_scenario[i,19] <- sum(dat.cf$T.A0>4)
descriptives_scenario[i,20] <- sum(dat.cf$T.A1>4)

colnames(descriptives_scenario)[17:20] <- c("n_at_risk_val0", "n_at_risk_val1", "n_at_risk_dat_cf0", "n_at_risk_dat_cf1")

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






