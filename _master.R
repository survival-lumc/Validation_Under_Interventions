#------------------------
#------------------------
#Preliminaries
#------------------------
#------------------------

#----
#set work directory


#----
#import libraries
library(tidyverse)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(xtable)
library(janitor)
library(survival)
library(timereg)
library(dynpred)
library(timeROC)

#----
#Import separate functions 
source("functions_discrimination.R")
source("functions_brier.R")

#------------------------
#------------------------
#Run the simulation in which data are generated and analysed using an Aalen additive hazards model 
#------------------------
#------------------------

#number of simulations
Nsim <- 1000

#sample size 
n=3000

#specify the scenario (+ plotting limits)
#scenario 1: base scenario: same dgm for development and validation data 
#scenario 2: different baseline hazard in development data: alpha.0 = - 1
#scenario 3: as scenario 1 plus noise added to L (L*) when making predictions in the validation set 
scenario<-1; censoring <- FALSE; cindex_ylim_low <- .48; cindex_ylim_high <- .625; auc_ylim_low <- 0.46; auc_ylim_high <- 0.70; calib_lim_low <- 0.40; calib_lim_high <- 0.95; calib_risk_lim_high <- 0.85; model <- " additive hazards model"
scenario<-2; censoring <- FALSE; cindex_ylim_low <- .48; cindex_ylim_high <- .625; auc_ylim_low <- 0.46; auc_ylim_high <- 0.70; calib_lim_low <- 0.45; calib_lim_high <- 0.95; calib_risk_lim_high <- 1; model <- " additive hazards model"
scenario<-3; censoring <- FALSE; cindex_ylim_low <- .48; cindex_ylim_high <- .625; auc_ylim_low <- 0.46; auc_ylim_high <- 0.70; calib_lim_low <- 0.40; calib_lim_high <- 0.95; calib_risk_lim_high <- 0.8; model <- " additive hazards model"

#perform analyses for Nsim simulated data sets of n individuals
source("simulation_addhaz.R")

#create results tables
source("results_tables.R")

#create results plots
source("results_plots.R")

#---
#print main results
descr_KM_plot
main_plots_together
appendix_plots_together
table_OErisk_wide
table_cindex_wide
table_auct_wide
table_ipa_wide

colMeans(descriptives_scenario)

#---
#save main plots in .png format

ggsave(descr_KM_plot,file=paste0("results_addhaz\\add_haz_marginal_risk_distribution","scenario_", scenario,".png"),width = 20/2, units = "cm")
ggsave(main_plots_together,file=paste0("results_addhaz\\add_haz_main_plots_together_","scenario_", scenario,".png"),width = 13, height = 20, units = "cm")
ggsave(appendix_plots_together,file=paste0("results_addhaz\\add_haz_appendix_plots_together_","scenario_", scenario,".png"),width = 13, height = 40/3, units = "cm")

#---
#save main tables

save(table_OErisk_wide,file=paste0("results_addhaz/table_OErisk_wide_","scenario_",scenario,".rda"))
save(table_cindex_wide,file=paste0("results_addhaz/table_cindex_wide_","scenario_",scenario,".rda"))
save(table_auct_wide,file=paste0("results_addhaz/table_auct_wide_","scenario_",scenario,".rda"))
save(table_ipa_wide,file=paste0("results_addhaz/table_ipa_wide_","scenario_",scenario,".rda"))

#---
#save raw result data
 
save(disc_cindex,file=paste0("results_addhaz/disc_cindex_","scenario_",scenario,".rda"))
save(disc_auct,file=paste0("results_addhaz/disc_auct_","scenario_",scenario,".rda"))
save(brier_raw,file=paste0("results_addhaz/brier_raw_","scenario_",scenario,".rda"))
save(brier_ipa,file=paste0("results_addhaz/brier_ipa_","scenario_",scenario,".rda"))
save(calib_risk0,file=paste0("results_addhaz/calib_risk0_","scenario_",scenario,".rda"))
save(calib_risk1,file=paste0("results_addhaz/calib_risk1_","scenario_",scenario,".rda"))
save(calib_risk0_group,file=paste0("results_addhaz/calib_risk0_group_","scenario_",scenario,".rda"))
save(calib_risk1_group,file=paste0("results_addhaz/calib_risk1_group_","scenario_",scenario,".rda"))
save(descriptives_scenario,file=paste0("results_addhaz/descriptives_","scenario_",scenario,".rda"))
save(descr_KM_curves_long,file=paste0("results_addhaz/descriptives_","scenario_",scenario,".rda"))

#---
#save separate plots in png

ggsave(cindex_plot_0,file=paste0("results_addhaz\\cindex_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(cindex_plot_1,file=paste0("results_addhaz\\cindex_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(auct_plot_0,file=paste0("results_addhaz\\auct_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(auct_plot_1,file=paste0("results_addhaz\\auct_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(brier_plot_0,file=paste0("results_addhaz\\brier_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(brier_plot_1,file=paste0("results_addhaz\\brier_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(ipa_plot_0,file=paste0("results_addhaz\\ipa_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(ipa_plot_1,file=paste0("results_addhaz\\ipa_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk_plot0,file=paste0("results_addhaz/calib_risk_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk_plot1,file=paste0("results_addhaz/calib_risk_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk0_group_plot,file=paste0("results_addhaz/calib_risk0_group_plot_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk1_group_plot,file=paste0("results_addhaz/calib_risk1_group_plot_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")

#---
#save plots in rda

save(descr_KM_plot,file=paste0("results_addhaz/descr_KM_plot_","scenario_",scenario,".rda"))
save(main_plots_together,file=paste0("results_addhaz\\main_plots_together_","scenario_", scenario,".rda"))
save(appendix_plots_together,file=paste0("results_addhaz\\appendix_plots_together_","scenario_", scenario,".rda"))

save(cindex_plot_0,file=paste0("results_addhaz\\cindex_plot0_","scenario_", scenario,".rda"))
save(cindex_plot_1,file=paste0("results_addhaz\\cindex_plot1_","scenario_", scenario,".rda"))
save(auct_plot_0,file=paste0("results_addhaz\\auct_plot0_","scenario_", scenario,".rda"))
save(auct_plot_1,file=paste0("results_addhaz\\auct_plot1_","scenario_", scenario,".rda"))
save(brier_plot_0,file=paste0("results_addhaz\\brier_plot0_","scenario_", scenario,".rda"))
save(brier_plot_1,file=paste0("results_addhaz\\brier_plot1_","scenario_", scenario,".rda"))
save(ipa_plot_0,file=paste0("results_addhaz\\ipa_plot0_","scenario_", scenario,".rda"))
save(ipa_plot_1,file=paste0("results_addhaz\\ipa_plot1_","scenario_", scenario,".rda"))
save(file=paste0("results_addhaz/calib_risk_plot0_","scenario_", scenario,".rda"),plot=calib_risk_plot0)
save(file=paste0("results_addhaz/calib_risk_plot1_","scenario_", scenario,".rda"),plot=calib_risk_plot1)
save(file=paste0("results_addhaz/calib_risk0_group_plot_","scenario_", scenario,".rda"),plot=calib_risk0_group_plot)
save(file=paste0("results_addhaz/calib_risk1_group_plot_","scenario_", scenario,".rda"),plot=calib_risk1_group_plot)

# #--------------
# # arrange descriptive KM plots from different scenarios together in one plot
# 
# descr_KM_plot1 <- descr_KM_plot  #after running scenario 1 (same as scenario 3)
# descr_KM_plot2 <- descr_KM_plot  #after running scenario 2
# KMplots_together <- ggarrange(descr_KM_plot1, descr_KM_plot2,
#                              nrow=1,
#                              common.legend=TRUE,legend="bottom")
# KMplots_together
# ggsave(KMplots_together,file=paste0("results_addhaz\\add_haz_KM_plots_together.png"),width = 20, units = "cm")
# 
# 
# #------------------------------
# #render latex code of main tables
# 
# xtable(table_OErisk_wide)
# xtable(table_cindex_wide)
# xtable(table_auct_wide)
# xtable(table_ipa_wide)

#------------------------
#------------------------
#Run the simulation in which data are generated and analysed using a Cox model 
#------------------------
#------------------------

#number of simulations
Nsim <- 1000

#sample size 
n=3000

# specify the scenario (+ plotting limits)
# scenario 1: base scenario same dgm for development and validation data (gamma.0=-1, gamma.L = 0.5, alpha.0 = -2, alpha.L = 0.5, alpha. A= -0.5, alpha.U= 0.5)
# scenario 2: different baseline hazard in development data (alpha.0 = - 1)
# scenario 3: as scenario 1 plus noise added to L (L*) when making predictions in the validation set 

scenario <- 1; cindex_ylim_low <- .48; cindex_ylim_high <- .675; auc_ylim_low <- 0.54; auc_ylim_high <- 0.75; calib_lim_low <- 0.1; calib_lim_high <- 0.9; calib_risk_lim_high <- 0.8 ;calib_events_lim_high <- (n/1000)*1300; model = " Cox model"
scenario <- 2; cindex_ylim_low <- .53; cindex_ylim_high <- .675; auc_ylim_low <- 0.54; auc_ylim_high <- 0.73; calib_lim_low <- 0.1; calib_lim_high <- 1.0; calib_risk_lim_high <- 1.0 ;calib_events_lim_high <- (n/1000)*3000; model = " Cox model"
scenario <- 3; cindex_ylim_low <- .52; cindex_ylim_high <- .675; auc_ylim_low <- 0.50; auc_ylim_high <- 0.70; calib_lim_low <- 0.1; calib_lim_high <- 0.9; calib_risk_lim_high <- 0.8 ;calib_events_lim_high <- (n/1000)*1300; model = " Cox model"

#Perform analyses for Nsim simulated data sets of n individuals
source("simulation_cox.R")

#create results tables
source("results_tables.R")

#create results plots
source("results_plots.R")

#---
#print main results
colMeans(descriptives_scenario)
descr_KM_plot
main_plots_together
appendix_plots_together
table_OErisk_wide
table_cindex_wide
table_auct_wide
table_ipa_wide

#---
#save main plots in .png format
#---
# Supplementary Figure A4 
ggsave(descr_KM_plot,file=paste0("results_cox\\cox_marginal_risk_distribution_","scenario_", scenario,".png"),width = 20/2, units = "cm")
# Supplementary Figure A5 (scenario 1), Figure A7 (scenario 2), Figure A9 (scenario 3)
ggsave(main_plots_together,file=paste0("results_cox\\cox_main_plots_together_","scenario_", scenario,".png"),width = 13, height = 20, units = "cm")
# Supplementary Figure A6 (scenario 1), Figure A8 (scenario 2), Figure A10 (scenario 3)
ggsave(appendix_plots_together,file=paste0("results_cox\\cox_appendix_plots_together_","scenario_", scenario,".png"),width = 13, height = 40/3, units = "cm")

#---
#save main tables

# Supplementary Table A2 (scenario 1), Table A3 (scenario 2), Table A4 (scenario 3)
save(table_OErisk_wide,file=paste0("results_cox/table_OErisk_wide_","scenario_",scenario,".rda"))
save(table_cindex_wide,file=paste0("results_cox/table_cindex_wide_","scenario_",scenario,".rda"))
save(table_auct_wide,file=paste0("results_cox/table_auct_wide_","scenario_",scenario,".rda"))
save(table_ipa_wide,file=paste0("results_cox/table_ipa_wide_","scenario_",scenario,".rda"))

#---
#save raw result data

save(disc_cindex,file=paste0("results_cox/disc_cindex_","scenario_",scenario,".rda"))
save(disc_auct,file=paste0("results_cox/disc_auct_","scenario_",scenario,".rda"))
save(brier_raw,file=paste0("results_cox/brier_raw_","scenario_",scenario,".rda"))
save(brier_ipa,file=paste0("results_cox/brier_ipa_","scenario_",scenario,".rda"))
save(calib_risk0,file=paste0("results_cox/calib_risk0_","scenario_",scenario,".rda"))
save(calib_risk1,file=paste0("results_cox/calib_risk1_","scenario_",scenario,".rda"))
save(calib_risk0_group,file=paste0("results_cox/calib_risk0_group_","scenario_",scenario,".rda"))
save(calib_risk1_group,file=paste0("results_cox/calib_risk1_group_","scenario_",scenario,".rda"))
save(descriptives_scenario,file=paste0("results_cox/descriptives_","scenario_",scenario,".rda"))
save(descr_KM_curves_long,file=paste0("results_cox/descriptives_","scenario_",scenario,".rda"))

#---
#save separate plots in png format

ggsave(cindex_plot_0,file=paste0("results_cox\\cindex_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(cindex_plot_1,file=paste0("results_cox\\cindex_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(auct_plot_0,file=paste0("results_cox\\auct_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(auct_plot_1,file=paste0("results_cox\\auct_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(brier_plot_0,file=paste0("results_cox\\brier_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(brier_plot_1,file=paste0("results_cox\\brier_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(ipa_plot_0,file=paste0("results_cox\\ipa_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(ipa_plot_1,file=paste0("results_cox\\ipa_plot1_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk_plot0, ,file=paste0("results_cox\\calib_risk_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk_plot1, ,file=paste0("results_cox\\calib_risk_plot0_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk0_group_plot, file=paste0("results_cox/calib_risk0_group_plot_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")
ggsave(calib_risk1_group_plot, file=paste0("results_cox/calib_risk1_group_plot_","scenario_", scenario,".png"),width = 20, height = 10, units = "cm")

#---
#save plots in rda format

 save(descr_KM_plot,file=paste0("results_cox/descr_KM_plot_","scenario_",scenario,".rda"))
 save(main_plots_together,file=paste0("results_cox\\main_plots_together_","scenario_", scenario,".rda"))
 save(appendix_plots_together,file=paste0("results_cox\\appendix_plots_together_","scenario_", scenario,".rda"))
 save(cindex_plot_0,file=paste0("results_cox\\cindex_plot0_","scenario_", scenario,".rda"))
 save(cindex_plot_1,file=paste0("results_cox\\cindex_plot1_","scenario_", scenario,".rda"))
 save(auct_plot_0,file=paste0("results_cox\\auct_plot0_","scenario_", scenario,".rda"))
 save(auct_plot_1,file=paste0("results_cox\\auct_plot1_","scenario_", scenario,".rda"))
 save(brier_plot_0,file=paste0("results_cox\\brier_plot0_","scenario_", scenario,".rda"))
 save(brier_plot_1,file=paste0("results_cox\\brier_plot1_","scenario_", scenario,".rda"))
 save(ipa_plot_0,file=paste0("results_cox\\ipa_plot0_","scenario_", scenario,".rda"))
 save(ipa_plot_1,file=paste0("results_cox\\ipa_plot1_","scenario_", scenario,".rda"))
 save(calib_risk_plot0, file=paste0("results_cox/calib_risk_plot0_","scenario_", scenario,".rda"))
 save(calib_risk_plot1, file=paste0("results_cox/calib_risk_plot1_","scenario_", scenario,".rda"))
 save(plot=calib_risk0_group_plot, file=paste0("results_cox/calib_risk0_group_plot_","scenario_", scenario,".rda"))
 save(plot=calib_risk1_group_plot, file=paste0("results_cox/calib_risk1_group_plot_","scenario_", scenario,".rda"))

#---
# arrange descriptive KM plots from different scenarios together in one plot

# descr_KM_plot1 <- descr_KM_plot  #after running scenario 1 (note same as scenario 3)
# descr_KM_plot2 <- descr_KM_plot  #after running scenario 2
# KMplots_together <- ggarrange(descr_KM_plot1, descr_KM_plot2,
#                              nrow=1,
#                              common.legend=TRUE,legend="bottom")
# KMplots_together
# ggsave(KMplots_together,file=paste0("results_cox\\cox_KM_plots_together.png"),width = 20, units = "cm")


#------------------------------
#render latex code of main tables

# xtable(table_OErisk_wide)
# xtable(table_cindex_wide)
# xtable(table_auct_wide)
# xtable(table_ipa_wide)
