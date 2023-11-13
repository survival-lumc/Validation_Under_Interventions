#---------------------
#---------------------
# Runs the simulations in which 
#     data are generated using a Cox model 
#     data are analysed using a Cox model (development of prediction model)
#     measures of predictive performance are calculated 
#     measures are stored in result matrices an arrays
#---------------------
#---------------------

#----
#seeds
set.seed(1)
seeds<-sample(1:1000000,Nsim)

#----
#storage for simulation results

#descriptives
descriptives_scenario <- matrix(nrow=Nsim, ncol=20)
colnames(descriptives_scenario) <- rep("NA",20)
descr_KM_curves_long <- matrix(nrow = 0, ncol = 4)
colnames(descr_KM_curves_long) = c("sim","time","type","risk")

#calibration
#[ , ,1]: predicted
#[ , ,2]: IPCW observed
#[ , ,3]: counterfactual true
#[ , ,4]: subset predicted
#[ , ,5]: subset observed
calib_risk0<-array(dim=c(Nsim,5,5))
calib_risk1<-array(dim=c(Nsim,5,5))
calib_risk0_group<-array(dim=c(Nsim,10,5))
calib_risk1_group<-array(dim=c(Nsim,10,5))

#discrimination
disc_cindex<-matrix(nrow=Nsim,ncol=10)
disc_auct<-matrix(nrow=Nsim,ncol=10)

#overall prediction error
brier_raw <- matrix(nrow=Nsim, ncol=8)
brier_ipa <- matrix(nrow=Nsim, ncol=8)

#----
#simulation loop

for(i in 1:Nsim){
  print(i)
  set.seed(seeds[i])
  
  #---
  #Simulate development and validation data and counterfactual validation data
  
  if(scenario==1|scenario=="5a"|scenario=="5b"|scenario=="6a"|scenario=="6b"|scenario=="6d"){
    source("dat_sim_cox_scenario1.R")
    dat.dev<-dat
    dat.long.dev<-dat.long
  
    source("dat_sim_cox_scenario1.R")
    dat.val<-dat
    dat.long.val<-dat.long
  
  }
  
  if(scenario==2){
    source("dat_sim_cox_scenario2.R")
    dat.dev<-dat
    dat.long.dev<-dat.long
    
    source("dat_sim_cox_scenario1.R")
    dat.val<-dat
    dat.long.val<-dat.long
    
  }
  
  if(scenario==3){
    source("dat_sim_cox_scenario1.R")
    dat.dev<-dat
    dat.long.dev<-dat.long
    
    source("dat_sim_cox_scenario3.R")
    dat.val<-dat
    dat.long.val<-dat.long

  }
  
  if(scenario=="4a"){
    source("dat_sim_cox_scenario1.R")
    dat.dev<-dat
    dat.long.dev<-dat.long
    
    source("dat_sim_cox_scenario4a.R")
    dat.val<-dat
    dat.long.val<-dat.long
    
  }
  
  if(scenario=="4b"){
    source("dat_sim_cox_scenario1.R")
    dat.dev<-dat
    dat.long.dev<-dat.long
    
    source("dat_sim_cox_scenario4b.R")
    dat.val<-dat
    dat.long.val<-dat.long
    
  }
  
  if(scenario=="6c"){
    source("dat_sim_cox_scenario1.R")
    dat.dev<-dat
    dat.long.dev<-dat.long
    
    source("dat_sim_cox_scenario6c.R")
    dat.val<-dat
    dat.long.val<-dat.long

  }
  
  #---
  #Simulate counterfactual validation data
  source("dat_sim_cox_counterfactual_scenario1.R")
  
  #---
  #Model development (uses dat.long.dev)
  source("analysis_cox_model_development.R")
  
  #---
  #Preparation for model validation: obtains quantities used for discrimination, calibration, Brier score (uses dat.val and dat.long.val)
  source("analysis_cox_model_validation.R")
  #source("analysis_cox_model_validation_weights_misspec.R")
  
  #---
  #Estimation of discrimination measures
  source("discrimination.R")
  
  #---
  #Estimation of calibration measures
  source("calibration.R")
  
  #Calculate descriptives for comparison of scenarios
  source("descriptives.R")
  
}
