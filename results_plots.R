###########################################################################################
###########################################################################################
# Creates the following plots summarizing simulation results:
#
# calib_risk0_group_plot and calib_risk1_group_plot: observed outcome proportions against mean estimated risks by time 5 within tenths of the estimated risks 
# cindex_plot_0 and cindex_plot_1: boxplots of cindex truncated at time 5 
# auct_plot_0 and auct_plot_1: boxplots of auct at time 5
# ipa_plot_0 and ipa_plot_1: boxplots of scaled brier at time 5
# main_plots_together:  combines the above 4x2 plots in one figure
#                       presented in the manuscript in Figures 3/4/5 (add hazard model scenario 1/2/3)
#                       presented in supplementary Figures A5/A7/A9 (Cox model scenario 1/2/3)
#
# calib_risk_plot0 and calib_risk_plot1: plots of estimated versus observed overall risks for time point 1,2,3,4,5
# brier_plot_0 and brier_plot_1: boxplots of (unscaled) brier at time 5
# appendix_plots_together:  combines the above 2x2 plots in one figure
#                           presented in supplementary Figures A1/A2/A3 (add hazard model scenario 1/2/3)
#                           presented in supplementary Figures A6/A8/A10 (Cox model scenario 1/2/3)
#
# descr_KM_plot: plot of the marginal risk distribution in development, validation, counterfactual never treated and always treated data
#                presented in the manuscript in Figure 2 (add hazards model)
#                presented in supplementary Figure A4 (Cox model)          
#
# The plots are saved in the master file
# 
###########################################################################################


#------------------------------
#------------------------------
#Summarize discrimination results: plots
#------------------------------
#------------------------------

# select estimators, colors and shapes of outliers to present in plots for cindex, auct, brier, ipa
selection_estimators <- c("true counterfactual","subset","IPCW")
selection_cols <- c("cf"="#000000","restr"="#D55E00","IPCW"="#009E73")

#---
#cindex

#restructuring of simulation results
disc_cindex0 <- disc_cindex[,c(1,9,2:4)]
disc_cindex1 <- disc_cindex[,c(5,10,6:8)]
colnames(disc_cindex0) <- c("true counterfactual", "subset", "Harrell's", "Uno's C", "IPCW")
colnames(disc_cindex1) <- c("true counterfactual", "subset", "Harrell's", "Uno's C", "IPCW")

sim=1:Nsim
dat_cindex0 <- cbind(disc_cindex0,sim)
dat_cindex1 <- cbind(disc_cindex1,sim)

disc_cindex_melt0<-melt(disc_cindex0,id=c("sim"), value.name = "cindex", varnames = c("sim","estimator"))
disc_cindex_melt1<-melt(disc_cindex1,id=c("sim"), value.name = "cindex", varnames = c("sim","estimator"))

cindex_plot_0<-ggplot(disc_cindex_melt0 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=cindex)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "C-index", x = ""
       #,title = "untreated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

cindex_plot_1<-ggplot(disc_cindex_melt1 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=cindex)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "C-index", x = ""
       #,title = "treated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

#---
#auct

#restructuring of simulation results
disc_auct0 <- disc_auct[,c(1,9,2:4)]
disc_auct1 <- disc_auct[,c(5,10,6:8)]
colnames(disc_auct0) <- c("true counterfactual", "subset", "unweighted", "KM weighted", "IPCW")
colnames(disc_auct1) <- c("true counterfactual", "subset", "unweighted", "KM weighted", "IPCW")

sim=1:Nsim
dat_auct0 <- cbind(disc_auct0,sim)
dat_auct1 <- cbind(disc_auct,sim)

disc_auct_melt0<-melt(disc_auct0,id=c("sim"), value.name = "auct", varnames = c("sim","estimator"))
disc_auct_melt1<-melt(disc_auct1,id=c("sim"), value.name = "auct", varnames = c("sim","estimator"))

auct_plot_0<-ggplot(disc_auct_melt0 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=auct)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "AUCt", x = ""
       #,title = "untreated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

auct_plot_1<-ggplot(disc_auct_melt1 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=auct)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "AUCt", x = ""
       #,title = "treated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

#---
#brier

brier0 <- brier_raw[,c(1,7,2,3)]
brier1 <- brier_raw[,c(4,8,5,6)]
colnames(brier0) <- c("true counterfactual", "subset", "KM weighted", "IPCW")
colnames(brier1) <- c("true counterfactual", "subset", "KM weighted", "IPCW")
sim=1:Nsim
dat_brier0 <- cbind(brier0,sim)
dat_brier1 <- cbind(brier1,sim)

brier_melt0<-melt(brier0,id=c("sim"), value.name = "brier", varnames = c("sim","estimator"))
brier_melt1<-melt(brier1,id=c("sim"), value.name = "brier", varnames = c("sim","estimator"))

brier_plot_0<-ggplot(brier_melt0 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=brier)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "Brier", x = ""
       #,title = "untreated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

brier_plot_1<-ggplot(brier_melt1 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=brier)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "Brier", x = ""
       #,title = "treated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

#---
#ipa

ipa0 <- brier_ipa[,c(1,7,2,3)]
ipa1 <- brier_ipa[,c(4,8,5,6)]
colnames(ipa0) <- c("true counterfactual", "subset", "KM weighted", "IPCW")
colnames(ipa1) <- c("true counterfactual", "subset", "KM weighted", "IPCW" )

sim=1:Nsim
dat_ipa0 <- cbind(ipa0,sim)
dat_ipa1 <- cbind(ipa,sim)

ipa_melt0<-melt(ipa0,id=c("sim"), value.name = "ipa", varnames = c("sim","estimator"))
ipa_melt1<-melt(ipa1,id=c("sim"), value.name = "ipa", varnames = c("sim","estimator"))

ipa_plot_0<-ggplot(ipa_melt0 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=ipa)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "Brier (%)", x = ""
       #,title = "untreated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

ipa_plot_1<-ggplot(ipa_melt1 %>% subset(., subset = estimator %in% selection_estimators), aes(x=estimator, y=ipa)) + 
  theme_bw()+
  geom_boxplot(col=selection_cols, outlier.shape = 20)+
  labs(y = "Brier (%)", x = ""
       #,title = "treated"
  )+
  theme(axis.ticks=element_blank(),
        #axis.text.x = element_text(angle = -45,vjust=-1,hjust=1.1)
        axis.text.x = element_blank() ) 

#------------------------------
#------------------------------
#Summarize calibration results: risk plot 
#------------------------------
#------------------------------

#reminder of content dim 3 'calib_risk0'  and 'calib_risk1'
#,,1 estimated
#,,2 observed / IPCW
#,,3 true / true cf
#,,4 estimated subset
#,,5 observed subset

#-----
#prepare data for plots

# estimated
lower_risk0_exp<-colMeans(calib_risk0[,,1])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,1])})
upper_risk0_exp<-colMeans(calib_risk0[,,1])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,1])})
# observed IPCW
lower_risk0_obs<-colMeans(calib_risk0[,,2])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,2])})
upper_risk0_obs<-colMeans(calib_risk0[,,2])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,2])})
# true cf
lower_risk0_true<-colMeans(calib_risk0[,,3])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,3])})
upper_risk0_true<-colMeans(calib_risk0[,,3])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,3])})
# subset expected
lower_risk0_exp_subset<-colMeans(calib_risk0[,,4])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,4])})
upper_risk0_exp_subset<-colMeans(calib_risk0[,,4])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,4])})
# subset observed
lower_risk0_obs_subset<-colMeans(calib_risk0[,,5])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,5])})
upper_risk0_obs_subset<-colMeans(calib_risk0[,,5])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk0[,x,5])})

calib_risk0_plot_data<-data.frame(time=1:5,timeA=(1:5)-0.2,timeB=(1:5)-0.1, timeC=(1:5)+0.1,timeD=(1:5)+0.2,
                                  est_risk0_exp=colMeans(calib_risk0[,,1]),lower_risk0_exp,upper_risk0_exp,
                                  est_risk0_obs=colMeans(calib_risk0[,,2]),lower_risk0_obs,upper_risk0_obs,
                                  est_risk0_true=colMeans(calib_risk0[,,3]),lower_risk0_true,upper_risk0_true,
                                  est_risk0_exp_subset=colMeans(calib_risk0[,,4]),lower_risk0_exp_subset,upper_risk0_exp_subset,
                                  est_risk0_obs_subset=colMeans(calib_risk0[,,5]),lower_risk0_obs_subset,upper_risk0_obs_subset
)

# estimated/expected
lower_risk1_exp<-colMeans(calib_risk1[,,1])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,1])})
upper_risk1_exp<-colMeans(calib_risk1[,,1])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,1])})
# observed IPCW
lower_risk1_obs<-colMeans(calib_risk1[,,2])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,2])})
upper_risk1_obs<-colMeans(calib_risk1[,,2])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,2])})
# true cf
lower_risk1_true<-colMeans(calib_risk1[,,3])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,3])})
upper_risk1_true<-colMeans(calib_risk1[,,3])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,3])})
# subset estimated/expected
lower_risk1_exp_subset<-colMeans(calib_risk1[,,4])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,4])})
upper_risk1_exp_subset<-colMeans(calib_risk1[,,4])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,4])})
# subset observed
lower_risk1_obs_subset<-colMeans(calib_risk1[,,5])-1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,5])})
upper_risk1_obs_subset<-colMeans(calib_risk1[,,5])+1.96*sapply(1:5,FUN=function(x){bias.mc(calib_risk1[,x,5])})

calib_risk1_plot_data<-data.frame(time=1:5,timeA=(1:5)-0.2,timeB=(1:5)-0.1, timeC=(1:5)+0.1,timeD=(1:5)+0.2,
                                  est_risk1_exp=colMeans(calib_risk1[,,1]),lower_risk1_exp,upper_risk1_exp,
                                  est_risk1_obs=colMeans(calib_risk1[,,2]),lower_risk1_obs,upper_risk1_obs,
                                  est_risk1_true=colMeans(calib_risk1[,,3]),lower_risk1_true,upper_risk1_true,
                                  est_risk1_exp_subset=colMeans(calib_risk1[,,4]),lower_risk1_exp,upper_risk1_exp_subset,
                                  est_risk1_obs_subset=colMeans(calib_risk1[,,5]),lower_risk1_obs,upper_risk1_obs_subset
)

#-----
#plot set-up functions

addpointtoplot <- function(dataset, varx, vary, vcol, vshape) { 
  list(
    geom_point(data=dataset, aes_string(x=varx, y=vary,colour=vcol, shape=vshape),size=1.5) 
  )
}

addsegmenttoplot <- function(dataset, varx, vary, varxend, varyend, vcol) { 
  list(
    geom_segment(data=dataset, aes_string(x=varx, y=vary, xend=varxend, yend=varyend,colour=vcol))
  )
}

#-----
# risk plot treated and untreated separate

cols   = c("true"="#000000","est"="grey", "obs"="#009E73", "est_subset"= "grey", "obs_subset"= "#D55E00")
shapes = c("true"=17,       "est"=17,     "obs"=17,        "est_subset"=15,         "obs_subset"=15)

calib_risk_plot0<-ggplot(calib_risk0_plot_data,aes(x=time,y=est_risk0_true))+theme_bw()+
  ylab("Risk")+xlab("Time")+
  addpointtoplot(calib_risk0_plot_data,"time","est_risk0_true",vcol='"true"',vshape = '"true"')+
  addsegmenttoplot(calib_risk0_plot_data,varx="time",vary="lower_risk0_true",varxend="time",varyend="upper_risk0_true",vcol='"true"')+
  addpointtoplot(calib_risk0_plot_data,"timeA","est_risk0_exp",vcol='"est"',vshape='"est"')+
  addsegmenttoplot(calib_risk0_plot_data,varx="timeA",vary="lower_risk0_exp",varxend="timeA",varyend="upper_risk0_exp",vcol='"est"')+
  addpointtoplot(calib_risk0_plot_data,"timeB","est_risk0_obs",vcol='"obs"',vshape='"obs"')+
  addsegmenttoplot(calib_risk0_plot_data,varx="timeB",vary="lower_risk0_obs",varxend="timeB",varyend="upper_risk0_obs",vcol='"obs"')+
  addpointtoplot(calib_risk0_plot_data,"timeC","est_risk0_exp_subset",vcol='"est"',vshape='"est_subset"')+
  addsegmenttoplot(calib_risk0_plot_data,varx="timeC",vary="lower_risk0_exp_subset",varxend="timeC",varyend="upper_risk0_exp_subset",vcol='"est"')+
  addpointtoplot(calib_risk0_plot_data,"timeD","est_risk0_obs_subset",vcol='"obs_subset"',vshape='"obs_subset"')+
  addsegmenttoplot(calib_risk0_plot_data,varx="timeD",vary="lower_risk0_obs_subset",varxend="timeD",varyend="upper_risk0_obs_subset",vcol='"obs_subset"')+
        #  
  scale_x_continuous(breaks=seq(1,5,1),limits=c(0.8,5.3))+
  scale_y_continuous(breaks=seq(0,calib_risk_lim_high,0.1),limits=c(0,calib_risk_lim_high))+
  scale_colour_manual(NULL,values=cols,labels=c(est="predicted", obs="estimated IPCW", true="true counterfactual",est_subset="predicted subset", obs_subset="estimated subset" ),
                      breaks=c("est","obs", "true","est_subset","obs_subset"))+
  scale_shape_manual(NULL,values=shapes,labels=c(est="predicted", obs="estimated IPCW", true="true counterfactual",est_subset="predicted subset", obs_subset="estimated subset" ),
                     breaks=c("est","obs", "true","est_subset","obs_subset"))+
  ggtitle("never treated")+
  theme(plot.title = element_text(hjust = 0.5))

calib_risk_plot1<-ggplot(calib_risk1_plot_data,aes(x=time,y=est_risk1_true))+theme_bw()+
  ylab("Risk")+xlab("Time")+
  addpointtoplot(calib_risk1_plot_data,"time","est_risk1_true",vcol='"true"',vshape = '"true"')+
  addsegmenttoplot(calib_risk1_plot_data,varx="time",vary="lower_risk1_true",varxend="time",varyend="upper_risk1_true",vcol='"true"')+
  addpointtoplot(calib_risk1_plot_data,"timeA","est_risk1_exp",vcol='"est"',vshape='"est"')+
  addsegmenttoplot(calib_risk1_plot_data,varx="timeA",vary="lower_risk1_exp",varxend="timeA",varyend="upper_risk1_exp",vcol='"est"')+
  addpointtoplot(calib_risk1_plot_data,"timeB","est_risk1_obs",vcol='"obs"',vshape='"obs"')+
  addsegmenttoplot(calib_risk1_plot_data,varx="timeB",vary="lower_risk1_obs",varxend="timeB",varyend="upper_risk1_obs",vcol='"obs"')+
  addpointtoplot(calib_risk1_plot_data,"timeC","est_risk1_exp_subset",vcol='"est"',vshape='"est_subset"')+
  addsegmenttoplot(calib_risk1_plot_data,varx="timeC",vary="lower_risk1_exp_subset",varxend="timeC",varyend="upper_risk1_exp_subset",vcol='"est"')+
  addpointtoplot(calib_risk1_plot_data,"timeD","est_risk1_obs_subset",vcol='"obs_subset"',vshape='"obs_subset"')+
  addsegmenttoplot(calib_risk1_plot_data,varx="timeD",vary="lower_risk1_obs_subset",varxend="timeD",varyend="upper_risk1_obs_subset",vcol='"obs_subset"')+
  #  
  scale_x_continuous(breaks=seq(1,5,1),limits=c(0.8,5.2))+
  scale_y_continuous(breaks=seq(0,calib_risk_lim_high,0.1),limits=c(0,calib_risk_lim_high))+
  scale_colour_manual(NULL,values=cols,labels=c(est="predicted", obs="estimated IPCW", true="true counterfactual",est_subset="predicted subset", obs_subset="estimated subset" ),
                      breaks=c("est","obs", "true","est_subset","obs_subset"))+
  scale_shape_manual(NULL,values=shapes,labels=c(est="predicted", obs="estimated IPCW", true="true counterfactual",est_subset="predicted subset", obs_subset="estimated subset" ),
                     breaks=c("est","obs", "true","est_subset","obs_subset"))+
  ggtitle("always treated")+
  theme(plot.title = element_text(hjust = 0.5))


#------------------------------
#------------------------------
#Summarize calibration results: Risks at time 5 in deciles of predicted risk, approach 1
#------------------------------
#------------------------------

#-----
#data for plots

calib_risk0_group_plot_data<-data.frame(est_risk0_group_exp=colMeans(calib_risk0_group[,,1],na.rm=T),
                                        est_risk0_group_obs=colMeans(calib_risk0_group[,,2],na.rm=T),
                                        est_risk0_group_true=colMeans(calib_risk0_group[,,3],na.rm=T),
                                        est_risk0_group_exp_subset=colMeans(calib_risk0_group[,,4],na.rm=T),
                                        est_risk0_group_obs_subset=colMeans(calib_risk0_group[,,5],na.rm=T))


calib_risk1_group_plot_data<-data.frame(est_risk1_group_exp=colMeans(calib_risk1_group[,,1],na.rm=T),
                                        est_risk1_group_obs=colMeans(calib_risk1_group[,,2],na.rm=T),
                                        est_risk1_group_true=colMeans(calib_risk1_group[,,3],na.rm=T),
                                        est_risk1_group_exp_subset=colMeans(calib_risk1_group[,,4],na.rm=T),
                                        est_risk1_group_obs_subset=colMeans(calib_risk1_group[,,5],na.rm=T))


#-----
#plot set-up

addpointtoplot <- function(dataset, varx, vary, vcol, vshape) { 
  list(
    geom_point(data=dataset, aes_string(x=varx, y=vary,colour=vcol, shape=vshape),size=1.5) 
  )
}

addlinetoplot <- function(dataset, varx, vary, vcol) { 
  list(
    geom_line(data=dataset, aes_string(x=varx, y=vary,colour=vcol)) 
  )
}

cols   = c("true" = "#000000", "obs" = "#009E73", "obs_subset" = "#D55E00")
shapes = c("true" = 19,        "obs" = 17,        "obs_subset" = 15)

#-----
#plot

calib_risk0_group_plot<-ggplot(calib_risk0_group_plot_data,aes(x=est_risk0_group_est,y=est_risk0_group_true))+
  theme_bw()+
  ylab("observed")+xlab("predicted")+
  #
  geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed", size=1)+
  #
  addpointtoplot(calib_risk0_group_plot_data,"est_risk0_group_exp","est_risk0_group_true",vcol='"true"',vshape='"true"')+
  addpointtoplot(calib_risk0_group_plot_data,"est_risk0_group_exp","est_risk0_group_obs",vcol='"obs"', vshape='"obs"')+
  addpointtoplot(calib_risk0_group_plot_data,"est_risk0_group_exp_subset","est_risk0_group_obs_subset",vcol='"obs_subset"', vshape='"obs_subset"')+
  addlinetoplot(calib_risk0_group_plot_data,"est_risk0_group_exp","est_risk0_group_true",vcol='"true"')+
  addlinetoplot(calib_risk0_group_plot_data,"est_risk0_group_exp","est_risk0_group_obs",vcol='"obs"')+
  addlinetoplot(calib_risk0_group_plot_data,"est_risk0_group_exp_subset","est_risk0_group_obs_subset",vcol='"obs_subset"')+
  #
  scale_x_continuous(breaks=seq(calib_lim_low,calib_lim_high,0.1),limits=c(calib_lim_low,calib_lim_high))+
  scale_y_continuous(breaks=seq(calib_lim_low,calib_lim_high,0.1),limits=c(calib_lim_low,calib_lim_high))+
  #
  scale_colour_manual(NULL,values=cols,labels=c(true="true counterfactual", obs_subset="estimated subset", obs="estimated IPCW"),
                      breaks=c("true","obs_subset","obs"))+
  scale_linetype_manual(NULL,values=cols,labels=c(true="true counterfactual", obs_subset="estimated subset", obs="estimated IPCW"),
                      breaks=c("true","obs_subset","obs"))+
  scale_shape_manual(NULL,values=shapes,labels=c(true="true counterfactual", obs="estimated IPCW", obs_subset="estimated subset" ),
                      breaks=c("true","obs_subset","obs"))+
  ggtitle("never treated")+
  theme(plot.title = element_text(hjust = 0.5 ))

calib_risk1_group_plot<-ggplot(calib_risk1_group_plot_data,aes(x=est_risk1_group_est,y=est_risk1_group_true))+theme_bw()+
  ylab("observed")+xlab("predicted")+
  #
  geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed", size=1)+
  #
  addpointtoplot(calib_risk1_group_plot_data,"est_risk1_group_exp","est_risk1_group_true",vcol='"true"',vshape='"true"')+
  addpointtoplot(calib_risk1_group_plot_data,"est_risk1_group_exp","est_risk1_group_obs",vcol='"obs"',vshape='"obs"')+
  addpointtoplot(calib_risk1_group_plot_data,"est_risk1_group_exp_subset","est_risk1_group_obs_subset",vcol='"obs_subset"', vshape='"obs_subset"')+
  addlinetoplot(calib_risk1_group_plot_data,"est_risk1_group_exp","est_risk1_group_true",vcol='"true"')+
  addlinetoplot(calib_risk1_group_plot_data,"est_risk1_group_exp","est_risk1_group_obs",vcol='"obs"')+
  addlinetoplot(calib_risk1_group_plot_data,"est_risk1_group_exp_subset","est_risk1_group_obs_subset",vcol='"obs_subset"')+
  #
  scale_x_continuous(breaks=seq(calib_lim_low,calib_lim_high,0.1),limits=c(calib_lim_low,calib_lim_high))+
  scale_y_continuous(breaks=seq(calib_lim_low,calib_lim_high,0.1),limits=c(calib_lim_low,calib_lim_high))+
  #
  scale_colour_manual(NULL,values=cols,labels=c(true="true counterfactual", obs_subset="estimated subset", obs="estimated IPCW"),
                      breaks=c("true","obs_subset","obs"))+
  scale_linetype_manual(NULL,values=cols,labels=c(true="true counterfactual", obs_subset="estimated subset", obs="estimated IPCW"),
                        breaks=c("true","obs_subset", "obs"))+
  scale_shape_manual(NULL,values=shapes,labels=c(true="true counterfactual", obs="estimated IPCW", obs_subset="estimated subset" ),
                     breaks=c("true","obs_subset", "obs"))+
  ggtitle("always treated")+
  theme(plot.title = element_text(hjust = 0.5))#to center the title

#------------------------------
#------------------------------
# arrange results together
#------------------------------
#------------------------------


main_plots_together <- ggarrange(
  calib_risk0_group_plot + theme(text = element_text(size = 8), legend.text = element_text(size = 10))
  , 
  calib_risk1_group_plot + theme(text = element_text(size = 8), legend.text = element_text(size = 10)),
  cindex_plot_0 + coord_cartesian(ylim = c(cindex_ylim_low, cindex_ylim_high)) + theme(text = element_text(size = 8), legend.text = element_text(size = 10), plot.margin = margin(b = 0)),  # Bottom margin
  cindex_plot_1 + coord_cartesian(ylim = c(cindex_ylim_low, cindex_ylim_high)) + theme(text = element_text(size = 8), legend.text = element_text(size = 10), plot.margin = margin(b = 0)),  # Bottom margin
  auct_plot_0 + coord_cartesian(ylim = c(auc_ylim_low,auc_ylim_high)) + theme(text = element_text(size = 8), legend.text = element_text(size = 10), plot.margin = margin(t = 0,  # Top margin
b = 0)),  # Bottom margin
  auct_plot_1 + coord_cartesian(ylim = c(auc_ylim_low,auc_ylim_high)) + theme(text = element_text(size = 8), legend.text = element_text(size = 10), plot.margin = margin(t = 0,  # Top margin
                                                                                                                        b = 0)),  # Bottom margin
  ipa_plot_0 + theme(text = element_text(size = 8), legend.text = element_text(size = 10), plot.margin = margin(t = 0)),  # Top margin
  ipa_plot_1 + theme(text = element_text(size = 8), legend.text = element_text(size = 10), plot.margin = margin(t = 0)), # top margin
  heights = c(1.55,1,1,1), nrow = 4,ncol=2, align = "v",
  #labels="AUTO",
  common.legend=TRUE,legend="bottom")

main_plots_together <- annotate_figure(main_plots_together, top = text_grob(paste0("Scenario ", scenario, model), size = 12))

appendix_plots_together <- ggarrange(
  calib_risk_plot0 + theme(text = element_text(size = 8), legend.text = element_text(size = 10)) + guides(col = guide_legend(nrow=2,byrow=TRUE)),
  calib_risk_plot1 + theme(text = element_text(size = 8), legend.text = element_text(size = 10)), 
  brier_plot_0 + theme(text = element_text(size = 10), plot.margin = margin(t = 0), legend.text = element_text(size = 8)),  # Top margin
  brier_plot_1 + theme(text = element_text(size = 10), plot.margin = margin(t = 0), legend.text = element_text(size = 8)), # top margin
  heights = c(1.1,0.95), nrow = 2,ncol=2, align = "v",
  #labels="AUTO",
  common.legend=TRUE,legend="bottom")

appendix_plots_together <- annotate_figure(appendix_plots_together, top = text_grob(paste0("Appendix scenario ", scenario, model), size = 12))


#-----------------------------
# descriptive plots
#-----------------------------

# marginal risk distribution

descr_KM_curves_average <- descr_KM_curves_long %>% group_by(type,time) %>% summarise(average_risk = mean(risk) ) 
descr_KM_curves_average$type <- factor(descr_KM_curves_average$type, levels=c("cf_untreated", "validation", "development", "cf_treated"))

cols=c("cf_untreated"= "#000000", "development"="#D55E00","validation"="#009E73", "cf_treated"="#000000")
ltypes=c("cf_untreated" = "solid", "development" = "dotdash", "validation" = "dashed", "cf_treated" = "dotted")

descr_KM_plot <- ggplot(data = descr_KM_curves_average, aes(x=time, y=average_risk, group = type, colour = type, linetype = type) )+
  geom_line(size=c(rep(1,51),rep(1,51),rep(1,51),rep(1,51)))+
  labs(x = "time",
       y = "cumulative risk of event",
       title = paste0("Scenario ", ifelse((scenario==1| scenario==3),"1 + 3", "2")))+
  scale_colour_manual(name=NULL, values=cols,     labels=c(cf_treated="counterfactual data treated", cf_untreated="counterfactual data untreated", development = "development data", validation = "validation data"))+
  scale_linetype_manual(name=NULL, values=ltypes, labels=c(cf_treated="counterfactual data treated", cf_untreated="counterfactual data untreated", development = "development data", validation = "validation data"))+
  ylim(0,0.8)


# L values (currently not used)

summary(dat.long$A)
dat.long$A <- as.character(dat.long$A)
dat.long$linetype <- ifelse(dat.long$A==0, "untreated", "treated")

cols=c("untreated" = "blue", "treated" = "red")
Lplot <- ggplot(data = dat.long, aes(x = visit, y = L, group = id, colour = linetype))+
  geom_line()+
  ylim(-23,33)+
  labs(x = "time",
       y = "L",
       title = paste0("Scenario ", scenario))+
  scale_colour_manual(name=NULL, values=cols, labels=c(treated="treated", untreated="untreated"))
  
