###########################################################################################
###########################################################################################
# Creates the following tables summarizing simulation results (mean, bias and Monte Carlo standard error)
#
# table_OErisk_wide: for calibration, contains ratios of observed vs expected risks by time point 5
# table_cindex_wide: for cindex by time point 5 
# table_auct_wide: for AUCt at time point 5
# table_ipa_wide: for scaled brier/ipa by time point 5
# table_combined_wide:  combines the above four tables into one table 
#                       presented in the main manuscript in Tables 1/2/3 (scenario 1/2/3 add hazard model) 
#                       presented in the supplementary file in eTables 4/5/6 (scenario 1/2/3 Cox model) 
#
# The tables are saved in the master file

###########################################################################################
###########################################################################################

#------------------------------
#------------------------------
#Functions used to summarize simulation results
#------------------------------
#------------------------------

#sp functions used later to display numbers to 1 or 3 decimal places
sp1=function(x){sprintf("%0.1f",x)}
sp3=function(x){sprintf("%0.3f",x)}

#function to calculate bias
bias=function(x,true){mean(x-true)}

#function to calculate MC error for the bias
bias.mc=function(x){sqrt(sum((x-mean(x))^2)/(Nsim*(Nsim-1)))}

#------------------------------
#------------------------------
#Table for calibration results
#------------------------------
#------------------------------
#----
#Calibration O/E ratio based on risk up to time 5
#Reminder of content dim 3 'calib_risk0'  and 'calib_risk1'
#,,1 estimated
#,,2 observed counterfactual / IPCW
#,,3 true / true cf
#,,4 estimated subset
#,,5 observed subset

calib_risk0_ratio <- cbind(calib_risk0[,5,3]/calib_risk0[,5,1], calib_risk0[,5,5]/calib_risk0[,5,4], calib_risk0[,5,2]/calib_risk0[,5,1])
colnames(calib_risk0_ratio) <- c("true", "subset", "counterfactual")

calib_risk0_ratio_summary<-sp3(colMeans(calib_risk0_ratio))
calib_risk0_ratio_bias<-sapply(1:3,FUN=function(x){bias(calib_risk0_ratio[,x],calib_risk0_ratio[,1])})
calib_risk0_ratio_biasmc<-sapply(1:3,FUN=function(x){bias.mc(calib_risk0_ratio[,x])})

calib_risk1_ratio <- cbind(calib_risk1[,5,3]/calib_risk1[,5,1], calib_risk1[,5,5]/calib_risk1[,5,4], calib_risk1[,5,2]/calib_risk1[,5,1])
colnames(calib_risk1_ratio) <- c("true", "subset", "counterfactual")

calib_risk1_ratio_summary<-sp3(colMeans(calib_risk1_ratio))
calib_risk1_ratio_bias<-sapply(1:3,FUN=function(x){bias(calib_risk1_ratio[,x],calib_risk1_ratio[,1])})
calib_risk1_ratio_biasmc<-sapply(1:3,FUN=function(x){bias.mc(calib_risk1_ratio[,x])})

table_OErisk0_wide <- rbind(calib_risk0_ratio_summary, paste0(sp3(calib_risk0_ratio_bias)," (",sp3(calib_risk0_ratio_biasmc),")"))
colnames(table_OErisk0_wide) <- c("true", "subset", "counterfactual")

table_OErisk1_wide <- rbind(calib_risk1_ratio_summary, paste0(sp3(calib_risk1_ratio_bias)," (",sp3(calib_risk1_ratio_biasmc),")"))
colnames(table_OErisk1_wide) <- c("true", "subset", "counterfactual")

table_OErisk0_wide[2,1] <- ""
table_OErisk1_wide[2,1] <- ""

table_OErisk_wide <- cbind(table_OErisk0_wide, table_OErisk1_wide)
rownames(table_OErisk_wide) <- c("mean", "bias (SE)")

#------------------------------
#------------------------------
#Tables for discrimination results
#------------------------------
#------------------------------

#---
#Cindex

#------------------------------
#restructuring of simulation results
disc_cindex0 <- disc_cindex[,c(1,9,2:4)]
disc_cindex1 <- disc_cindex[,c(5,10,6:8)]
colnames(disc_cindex0) <- c("true", "subset", "unweighted", "KM weighted", "counterfactual") #note: Harrell's = unweighted, Uno's C = KM weighted
colnames(disc_cindex1) <- c("true", "subset", "unweighted", "KM weighted", "counterfactual")

#------------------------------
#Select which estimators you want to see the results of

selection_estimators <- c("true","subset","counterfactual")

sim=1:Nsim
dat_cindex0 <- cbind(disc_cindex0,sim)
dat_cindex1 <- cbind(disc_cindex1,sim)

disc_cindex_melt0<-melt(disc_cindex0,id=c("sim"), value.name = "cindex", measure.vars = c("ignore treatment", "subset", "Harrell's", "Uno's C", "counterfactual"), varnames = c("sim","estimator"))
truth_cindex0 <- rep(disc_cindex0[,1], length(colnames(disc_cindex0)))
disc_cindex_melt0<-cbind(disc_cindex_melt0,truth_cindex0) #append the true values

table_cindex0 <-  disc_cindex_melt0 %>% 
                  subset(., subset = estimator %in% selection_estimators) %>%
                  group_by(estimator) %>%
                  summarise(mean = sp3(mean(cindex)),
                                bias = sp3(bias(cindex, truth_cindex0)),
                                MCerror = sp3(bias.mc(cindex))) %>% 
                  t() %>%
                  as.data.frame()


table_cindex0_wide <- disc_cindex_melt0 %>% 
                             subset(., subset = estimator %in% selection_estimators) %>%
                             group_by(estimator) %>%
                             summarise(mean = sp3(mean(cindex)),
                             biasMC = paste0(sp3(bias(cindex, truth_cindex0))," (", sp3(bias.mc(cindex)),")" )) %>%
                             t() %>%
                             as.data.frame()
  
disc_cindex_melt1<-melt(disc_cindex1,id=c("sim"), value.name = "cindex", varnames = c("sim","estimator"))
truth_cindex1 <- rep(disc_cindex1[,1], length(colnames(disc_cindex1)))
disc_cindex_melt1<-cbind(disc_cindex_melt1,truth_cindex1)

table_cindex1 <-  disc_cindex_melt1 %>% 
                  subset(., subset = estimator %in% selection_estimators) %>%
                  group_by(estimator) %>%
                  summarise(mean = sp3(mean(cindex)),
                            bias = sp3(bias(cindex, truth_cindex1)),
                            MCerror = sp3(bias.mc(cindex))) %>% 
                  t() %>%
                  as.data.frame()

table_cindex1_wide <- disc_cindex_melt1 %>% 
                      subset(., subset = estimator %in% selection_estimators) %>%
                      group_by(estimator) %>%
                      summarise(mean = sp3(mean(cindex)),
                      biasMC = paste0(sp3(bias(cindex, truth_cindex1))," (", sp3(bias.mc(cindex)),")" )) %>%
                      t() %>%
                      as.data.frame()

  
#make prettier
table_cindex1[4,1] <- ""
table_cindex1[3,1] <- ""
table_cindex0[4,1] <- ""
table_cindex0[3,1] <- ""

table_cindex0_wide[3,1] <- ""
table_cindex1_wide[3,1] <- ""

#combine never treated and always treated in one table
table_cindex <- cbind(table_cindex0, table_cindex1)
colnames(table_cindex) <- table_cindex[1,]
table_cindex <- table_cindex[-1,]

#make table wide (ie bias and SE next to each other)
table_cindex_wide <- cbind(table_cindex0_wide, table_cindex1_wide)
colnames(table_cindex_wide) <- table_cindex_wide[1,]
table_cindex_wide <- table_cindex_wide[-1,]
rownames(table_cindex_wide) <- c("mean", "bias (SE)")

#---
#C/D AUCt

#------------------------------
#restructuring of simulation results
disc_auct0 <- disc_auct[,c(1,9,2:4)]
disc_auct1 <- disc_auct[,c(5,10,6:8)]
colnames(disc_auct0) <- c("true", "subset", "unweighted", "KM weighted", "counterfactual")
colnames(disc_auct1) <- c("true", "subset", "unweighted", "KM weighted", "counterfactual")

sim=1:Nsim
dat_auct0 <- cbind(disc_auct0,sim)
dat_auct1 <- cbind(disc_auct,sim)

disc_auct_melt0<-melt(disc_auct0,id=c("sim"), value.name = "auct", varnames = c("sim","estimator"))
truth_auct0 <- rep(disc_auct0[,1], length(colnames(disc_auct0)))
disc_auct_melt0<-cbind(disc_auct_melt0,truth_auct0) #append the true values

table_auct0 <-  disc_auct_melt0 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(auct)),
            bias = sp3(bias(auct, truth_auct0)),
            MCerror = sp3(bias.mc(auct))) %>% 
  t() %>%
  as.data.frame()

table_auct0_wide <- disc_auct_melt0 %>% 
                    subset(., subset = estimator %in% selection_estimators) %>%
                    group_by(estimator) %>%
                    summarise(mean = sp3(mean(auct)),
                              biasMC = paste0(sp3(bias(auct, truth_auct0))," (", sp3(bias.mc(auct)),")" )) %>%
                    t() %>%
                    as.data.frame()

disc_auct_melt1<-melt(disc_auct1,id=c("sim"), value.name = "auct", varnames = c("sim","estimator"))
truth_auct1 <- rep(disc_auct1[,1], length(colnames(disc_auct1)))
disc_auct_melt1<-cbind(disc_auct_melt1,truth_auct1) #append the true values

table_auct1 <-  disc_auct_melt1 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(auct)),
            bias = sp3(bias(auct, truth_auct1)),
            MCerror = sp3(bias.mc(auct))) %>% 
  t() %>%
  as.data.frame()

table_auct1_wide <- disc_auct_melt1 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(auct)),
            biasMC = paste0(sp3(bias(auct, truth_auct1))," (", sp3(bias.mc(auct)),")" )) %>%
  t() %>%
  as.data.frame()

#make prettier
table_auct0[4,1] <- ""
table_auct0[3,1] <- ""
table_auct1[3,1] <- ""
table_auct1[4,1] <- ""

table_auct0_wide[3,1] <- ""
table_auct1_wide[3,1] <- ""

#combine never treated and always treated in one table
table_auct <- cbind(table_auct0, table_auct1)
colnames(table_auct) <- table_auct[1,]
table_auct <- table_auct[-1,]

#make table wide (ie bias and SE next to each other)
table_auct_wide <- cbind(table_auct0_wide, table_auct1_wide)
colnames(table_auct_wide) <- table_auct_wide[1,]
table_auct_wide <- table_auct_wide[-1,]
rownames(table_auct_wide) <- c("mean", "bias (SE)")


#------------------------------
#------------------------------
#Tables for Brier results
#------------------------------
#------------------------------

#---
#Brier

brier0 <- brier_raw[,c(1,7,2,3)]
brier1 <- brier_raw[,c(4,8,5,6)]
colnames(brier0) <- c("true", "subset", "KM weighted", "counterfactual")
colnames(brier1) <- c("true", "subset", "KM weighted", "counterfactual")
sim=1:Nsim
dat_brier0 <- cbind(brier0,sim)
brier_melt0<-melt(brier0,id=c("sim"), value.name = "brier", varnames = c("sim","estimator"))
truth_brier0 <- rep(brier0[,1], length(colnames(brier0)))
brier_melt0<-cbind(brier_melt0,truth_brier0) #append the true values

dat_brier1 <- cbind(brier1,sim)
brier_melt1<-melt(brier1,id=c("sim"), value.name = "brier", varnames = c("sim","estimator"))
truth_brier1 <- rep(brier1[,1], length(colnames(brier1)))
brier_melt1<-cbind(brier_melt1,truth_brier1) #append the true values

table_brier0 <-  brier_melt0 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(brier)),
            bias = sp3(bias(brier, truth_brier0)),
            MCerror = sp3(bias.mc(brier))) %>% 
  t() %>%
  as.data.frame()

table_brier0_wide <- brier_melt0 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(brier)),
            biasMC = paste0(sp3(bias(brier, truth_brier0))," (", sp3(bias.mc(brier)),")" )) %>%
  t() %>%
  as.data.frame()


table_brier1 <-  brier_melt1 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(brier)),
            bias = sp3(bias(brier, truth_brier1)),
            MCerror = sp3(bias.mc(brier))) %>% 
  t() %>%
  as.data.frame()

table_brier1_wide <- brier_melt1 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(brier)),
            biasMC = paste0(sp3(bias(brier, truth_brier1))," (", sp3(bias.mc(brier)),")" )) %>%
  t() %>%
  as.data.frame()

#make prettier
table_brier0[4,1] <- ""
table_brier0[3,1] <- ""
table_brier1[3,1] <- ""
table_brier1[4,1] <- ""

table_brier0_wide[3,1] <- ""
table_brier1_wide[3,1] <- ""

#combine never treated and always treated into one table
table_brier <- cbind(table_brier0, table_brier1)
colnames(table_brier) <- table_brier[1,]
table_brier <- table_brier[-1,]

#make table wide (ie bias and SE next to each other)
table_brier_wide <- cbind(table_brier0_wide, table_brier1_wide)
colnames(table_brier_wide) <- table_brier_wide[1,]
table_brier_wide <- table_brier_wide[-1,]
rownames(table_brier_wide) <- c("mean", "bias (SE)")

#---
#IPA (index of prediction accuracy or scaled Brier)

ipa0 <- brier_ipa[,c(1,7,2,3)]
ipa1 <- brier_ipa[,c(4,8,5,6)]
colnames(ipa0) <- c("true", "subset", "KM weighted", "counterfactual")
colnames(ipa1) <- c("true", "subset", "KM weighted", "counterfactual")
sim=1:Nsim
dat_ipa0 <- cbind(ipa0,sim)
ipa_melt0<-melt(ipa0,id=c("sim"), value.name = "ipa", varnames = c("sim","estimator"))
truth_ipa0 <- rep(ipa0[,1], length(colnames(ipa0)))
ipa_melt0<-cbind(ipa_melt0,truth_ipa0) #append the true values

dat_ipa1 <- cbind(ipa1,sim)
ipa_melt1<-melt(ipa1,id=c("sim"), value.name = "ipa", varnames = c("sim","estimator"))
truth_ipa1 <- rep(ipa1[,1], length(colnames(ipa1)))
ipa_melt1<-cbind(ipa_melt1,truth_ipa1) #append the true values

table_ipa0 <-  ipa_melt0 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(ipa)),
            bias = sp3(bias(ipa, truth_ipa0)),
            MCerror = sp3(bias.mc(ipa))) %>% 
  t() %>%
  as.data.frame()

table_ipa0_wide <- ipa_melt0 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(ipa)),
            biasMC = paste0(sp3(bias(ipa, truth_ipa0))," (", sp3(bias.mc(ipa)),")" )) %>%
  t() %>%
  as.data.frame()

table_ipa1 <-  ipa_melt1 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(ipa)),
            bias = sp3(bias(ipa, truth_ipa1)),
            MCerror = sp3(bias.mc(ipa))) %>% 
  t() %>%
  as.data.frame()

table_ipa1_wide <- ipa_melt1 %>% 
  subset(., subset = estimator %in% selection_estimators) %>%
  group_by(estimator) %>%
  summarise(mean = sp3(mean(ipa)),
            biasMC = paste0(sp3(bias(ipa, truth_ipa1))," (", sp3(bias.mc(ipa)),")" )) %>%
  t() %>%
  as.data.frame()

#make prettier
table_ipa0[4,1] <- ""
table_ipa0[3,1] <- ""
table_ipa1[3,1] <- ""
table_ipa1[4,1] <- ""

table_ipa0_wide[3,1] <- ""
table_ipa1_wide[3,1] <- ""

#combine never treated and always treated into one table
table_ipa <- cbind(table_ipa0, table_ipa1)
colnames(table_ipa) <- table_ipa[1,]
table_ipa <- table_ipa[-1,]

#wide format with bias and MC error next to each other
table_ipa_wide <- cbind(table_ipa0_wide, table_ipa1_wide)
colnames(table_ipa_wide) <- table_ipa_wide[1,]
table_ipa_wide <- table_ipa_wide[-1,]
rownames(table_ipa_wide) <- c("mean", "bias (SE)")

#------------------------------
#------------------------------
#Combine the above tables on OE ratio, cindex, auct and ipa into one table
#------------------------------
#------------------------------

table_combined_wide <- rbind(table_OErisk_wide, table_cindex_wide, table_auct_wide, table_ipa_wide)


 