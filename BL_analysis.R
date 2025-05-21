# Baseline differences: functionally stable group (FSG) vs functional loss group (FLG)
# Max Korbmacher, 03.12.2024
#
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# define data path
datapath = "/Users/max/Documents/Local/MS/results/"
# read data
data = read.csv(paste(datapath,"interrim_data.csv",sep=""))
# pkgs
library(dplyr)
library(datawizard)
library(rlist)
library(WebPower)
library(naniar)
library(ggplot2)
#
data = data_unique(data = data, select = eid)
#
# missingness
missi = data
names(missi) = gsub("Age_BL_OFAMS", "Age", names(missi))
names(missi) = gsub("baselineV", "Lesion volume", names(missi))
names(missi) = gsub("baselineC", "Lesion count", names(missi))
names(missi) = gsub("TotalVol", "Brain volume", names(missi))
names(missi) = gsub("TIV", "Intracranial volume", names(missi))
names(missi) = gsub("geno", "HLA-DRB1 carrier", names(missi))
names(missi) = gsub("CH3L.1..mg.ml..mean", "Chitinase-3 like-protein-1 mg/ml", names(missi))
names(missi) = gsub("BL_BMI", "Body Mass Index", names(missi))
names(missi) = gsub("edss_baseline", "EDSS", names(missi))
names(missi) = gsub("Vit_A_0", "Vitamin A", names(missi))
names(missi) = gsub("Vit_D_0", "Vitamin D", names(missi))
names(missi) = gsub("Vit_E_0", "Vitamin ", names(missi))
names(missi) = gsub("PF", "Physical functioning", names(missi))
names(missi) = gsub("RF", "Role-physical", names(missi))
names(missi) = gsub("BP", "Bodily pain", names(missi))
names(missi) = gsub("GH", "General health", names(missi))
names(missi) = gsub("VT", "Vitality", names(missi))
names(missi) = gsub("SF", "Social functioning", names(missi))
names(missi) = gsub("RE", "Role-emotional", names(missi))
names(missi) = gsub("MH", "Mental health", names(missi))
names(missi) = gsub("BAG_c", "Brain age gap", names(missi))
names(missi) = gsub("NfL..pg.ml.", "Neurofilament light chain", names(missi))
names(missi) = gsub("BAG_c", "Brain age gap", names(missi))
names(missi) = gsub("relapses_12mnths_before_baseline", "Relapses at baseline", names(missi)) # 12 months prior
names(missi) = gsub("smoking_OFAMS", "Smoking status", names(missi))
names(missi) = gsub("sex", "Sex", names(missi))
names(missi) = gsub("Treatment_OFAMS", "Omega-3 Intervention", names(missi))
names(missi) = gsub("age","Age",names(missi))
names(missi) = gsub("edss","Expanded Disability Status Scale (EDSS)",names(missi))
names(missi) = gsub("PASAT","Paced Auditory Serial Addition Test (PASAT)",names(missi))
names(missi) = gsub("Omega3_suppl","Omega 3 supplement received",names(missi))
missi=missi[order(names(missi)),] # order by name
missi = gg_miss_var(missi%>%dplyr::select(-c(eid,X,CLG,FLG)))
ggsave(paste(datapath,"missingness.pdf",sep=""),missi)
rm(missi)
#
#
################ 1. FUNCTIONAL DECLINE ####
# 1.1 POWER ####
# power is calculated from another script (power.R or power_optimized.R for parallelisation)
# we use some descriptives on power and model fit
model_info = read.csv(paste(datapath,"EDSS_model_info_optimised.csv",sep=""))
model_info = model_info %>% dplyr::filter(Formula != "~") %>% dplyr::filter(Formula != "FLG")
pwr = read.csv(paste(datapath,"power.csv",sep=""))
#pwr$Formula = gsub("FLG ~ ","",pwr$Formula)
pwr = merge(model_info,pwr,by="Formula")
write.csv(x = model_info, paste(datapath,"model_info.csv",sep=""))
rm(model_info)
# all models
pwr$Power = pwr$pwr
hist(pwr$Power)
hist(pwr$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of all models ran
sum(is.na(pwr$Power)) # models with convergence errors disallowing power calculation
sum(ifelse(pwr$McFaddenR2 > 0.1, 1, 0))/nrow(pwr)
#
# high enough power models
pw = pwr
pwr = pwr %>% filter(Power >= .8)
hist(pwr$Power)
hist(pwr$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of  models ran
nrow(pwr) / nrow(pw) # proportion of n
#
# high enough power and excellent model fit
pwr = pwr %>% filter(McFaddenR2 >= .2)
hist(pwr$Power)
hist(pwr$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of all models ran
nrow(pwr) / nrow(pw) # proportion of n
#
# excellent model fit
pwr = pw %>% filter(McFaddenR2 >= .2)
hist(pwr$Power)
hist(pwr$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of  models ran
nrow(pwr) / nrow(pw) # proportion of n
rm(pwr)
#
#
#
#
# 1.2 DEMOGRAPHICS: BL differences between FLG and non-FLG patients ####
# Differences in continuous scores ####
demotab=function(var){
  M_FLG = mean(unlist(na.omit(data%>%filter(FLG == 1))[var]))
  SD_FLG = sd(unlist(na.omit(data%>%filter(FLG == 1))[var]))
  M_FSG = mean(unlist(na.omit(data%>%filter(FLG == 0))[var]))
  SD_FSG = sd(unlist(na.omit(data%>%filter(FLG == 0))[var]))
  FLG.imp = lm(formula(paste(var,"~ FLG + sex + age")), data)
  B = summary(FLG.imp)$coefficients[2]
  SE = summary(FLG.imp)$coefficients[2,2]
  CI_low = B-1.96*SE
  CI_high = B+1.96*SE
  t = summary(FLG.imp)$coefficients[2,3]
  p = summary(FLG.imp)$coefficients[2,4]
  return(data.frame(var,M_FLG,SD_FLG,M_FSG,SD_FSG,B,SE,CI_low, CI_high, t,p))
}
contvar = c("DISEASE_DURATION","age_gap","BAG_c","baselineC","baselineV","edss","PASAT",
            "relapses_12mnths_before_baseline","BL_BMI",
            "NfL..pg.ml.",
            "Vit_A_0", "Vit_D_0","Vit_E_0", "BP", "GH", "MH",
            "PF", "SF", "VT") # constants: "RE", "RF"
demg = list()
for (i in 1:length(contvar)){
  demg[[i]] = demotab(contvar[i])
}

demg = list.rbind(demg)
# apply the whole spiel for age as well 
demotab2=function(var){
  M_FLG = mean(unlist(na.omit(data%>%filter(FLG == 1))[var]))
  SD_FLG = sd(unlist(na.omit(data%>%filter(FLG == 1))[var]))
  M_FSG = mean(unlist(na.omit(data%>%filter(FLG == 0))[var]))
  SD_FSG = sd(unlist(na.omit(data%>%filter(FLG == 0))[var]))
  FLG.imp = lm(formula(paste(var,"~ FLG + sex")), data)
  B = summary(FLG.imp)$coefficients[2]
  SE = summary(FLG.imp)$coefficients[2,2]
  CI_low = B-1.96*SE
  CI_high = B+1.96*SE
  t = summary(FLG.imp)$coefficients[2,3]
  p = summary(FLG.imp)$coefficients[2,4]
  return(data.frame(var,M_FLG,SD_FLG,M_FSG,SD_FSG,B,SE,CI_low, CI_high, t,p))
}

demg = rbind(demotab2("age"),(demg))
demg[2:(ncol(demg)-1)] = round(demg[2:(ncol(demg)-1)], 2)
demg$p = round(demg$p,4)
demg$FSG = paste(demg$M_FSG," (",demg$SD_FSG,")",sep="")
demg$FLG = paste(demg$M_FLG," (",demg$SD_FLG,")",sep="")
demg$Difference = paste(demg$B," (",demg$SE,")",sep="")
demg$CI = paste(round(demg$CI_low,2), "; ",round(demg$CI_high,2),sep="")
demg = demg %>% dplyr::select(FLG,FSG,Difference,CI,t,p)
write.csv(x = demg, paste(datapath,"descriptive_continuous.csv",sep=""))
#
# check also median and mad for EDSS
median(unlist(na.omit(data%>%filter(FLG == 1))["edss"]))
mad(unlist(na.omit(data%>%filter(FLG == 1))["edss"]))
median(unlist(na.omit(data%>%filter(FLG == 0))["edss"]))
mad(unlist(na.omit(data%>%filter(FLG == 0))["edss"]))
#
#
# Differences in frequencies ####
# FLG 1 = FLG
# sex: 1=female, 0=male
binaries = c("sex", "smoking_OFAMS","Omega3_suppl", "geno","Current_DMT")
res=list()
for (i in 1:length(binaries)){
  # how the sausage is made, but obsolete for OR, CI, p:
  a = table(data.frame(data["FLG"],data[paste(binaries[i])]))[4] # both 1
  b = table(data.frame(data["FLG"],data[paste(binaries[i])]))[1] # both 0
  c = table(data.frame(data["FLG"],data[paste(binaries[i])]))[3] # only FLG 1
  d = table(data.frame(data["FLG"],data[paste(binaries[i])]))[2] # only var of interest 1
  # OR = (a * b) / (c * d) # resulting odds ratio
  # SE = sqrt(1/a + 1/b + 1/c + 1/d)
  # CI_low = exp(log(OR)-1.96*SE)
  # CI_high = exp(log(OR)+1.96*SE)
  # data.frame(OR,SE,CI_low, CI_high)
  res[[i]] = c(a,b,c,d,(fisher.test(table(data.frame(data["FLG"],data[paste(binaries[i])])))$estimate),
    unlist(fisher.test(table(data.frame(data["FLG"],data[paste(binaries[i])])))[2]),
    unlist(fisher.test(table(data.frame(data["FLG"],data[paste(binaries[i])])))[1]))
}
res = data.frame(list.rbind(res))
rownames(res) = binaries
res$OR=paste(round(res$odds.ratio,2)," [",round(res$conf.int1,2),",",round(res$conf.int2,2),"]",sep="")
res = res %>% dplyr::select(V1,V2,V3,V4,OR,p.value)
names(res) = c("FLG & Group membership = yes","FLG & Group membership = no", 
               "FSG & Group membership = yes","FSG & Group membership = no",
               "OR","p")
write.csv(x = res, paste(datapath,"descriptive_binary.csv",sep=""))
#
#
################ 2. COGNITIVE DECLINE ####
# 2.1 POWER ####
# power is calculated from another script (power.R or power_optimized.R for parallelisation)
# we use some descriptives on power and model fit
model_info = read.csv(paste(datapath,"PASAT_model_info_optimised.csv",sep=""))
model_info = model_info %>% dplyr::filter(Formula != "~") %>% dplyr::filter(Formula != "CLG")
pwr = read.csv(paste(datapath,"CLG_power.csv",sep=""))
#pwr$Formula = gsub("CLG ~ ","",pwr$Formula)
pwr = merge(model_info,pwr,by="Formula")
rm(model_info)
# all models
pwr$Power = pwr$pwr
hist(pwr$Power)
hist(pwr$McFaddenR2)
median(na.omit(pwr$Power))
mad(na.omit(pwr$Power))
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of all models ran
sum(is.na(pwr$Power)) # models with convergence errors disallowing power calculation
sum(ifelse(pwr$McFaddenR2 > 0.1, 1, 0))/nrow(pwr)
#
# high enough power models
pw = pwr
pwr %>% filter(Power >= .8)
hist(pw$Power)
hist(pw$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of  models ran
nrow(pwr) / nrow(pw) # proportion of n
#
# high enough power and excellent model fit
pwr = pwr %>% filter(McFaddenR2 >= .2)
hist(pwr$Power)
hist(pwr$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of all models ran
nrow(pwr) / nrow(pw) # proportion of n
#
# excellent model fit
pwr = pw %>% filter(McFaddenR2 >= .2)
hist(pwr$Power)
hist(pwr$McFaddenR2)
paste("Power =",round(median(na.omit(pwr$Power)),2),"±",round(mad(na.omit(pwr$Power)),2))
paste("pR2 =",round(median(pwr$McFaddenR2),2),"±",round(mad(pwr$McFaddenR2),2))
paste("AUC =",round(median(pwr$AUC),2),"±",round(mad(pwr$AUC),2))
nrow(pwr) # number of  models ran
nrow(pwr) / nrow(pw) # proportion of n
rm(pwr)
#
#
#
#
# 2.2 DEMOGRAPHICS: BL differences between CLG and non-CLG patients ####
# Differences in continuous scores ####

demotab=function(var){
  M_CLG = mean(unlist(na.omit(data%>%filter(CLG == 1))[var]))
  SD_CLG = sd(unlist(na.omit(data%>%filter(CLG == 1))[var]))
  M_CSIG = mean(unlist(na.omit(data%>%filter(CLG == 0))[var]))
  SD_CSIG = sd(unlist(na.omit(data%>%filter(CLG == 0))[var]))
  CLG.imp = lm(formula(paste(var,"~ CLG + sex + age")), data)
  B = summary(CLG.imp)$coefficients[2]
  SE = summary(CLG.imp)$coefficients[2,2]
  CI_low = B-1.96*SE
  CI_high = B+1.96*SE
  t = summary(CLG.imp)$coefficients[2,3]
  p = summary(CLG.imp)$coefficients[2,4]
  return(data.frame(var,M_CLG,SD_CLG,M_CSIG,SD_CSIG,B,SE,CI_low, CI_high, t,p))
}

demg = list()
for (i in 1:length(contvar)){
  demg[[i]] = demotab(contvar[i])
}
demg = list.rbind(demg)
# apply the whole spiel for age as well 
demotab2=function(var){
  M_CLG = mean(unlist(na.omit(data%>%filter(CLG == 1))[var]))
  SD_CLG = sd(unlist(na.omit(data%>%filter(CLG == 1))[var]))
  M_CSIG = mean(unlist(na.omit(data%>%filter(CLG == 0))[var]))
  SD_CSIG = sd(unlist(na.omit(data%>%filter(CLG == 0))[var]))
  CLG.imp = lm(formula(paste(var,"~ CLG + sex")), data)
  B = summary(CLG.imp)$coefficients[2]
  SE = summary(CLG.imp)$coefficients[2,2]
  CI_low = B-1.96*SE
  CI_high = B+1.96*SE
  t = summary(CLG.imp)$coefficients[2,3]
  p = summary(CLG.imp)$coefficients[2,4]
  return(data.frame(var,M_CLG,SD_CLG,M_CSIG,SD_CSIG,B,SE,CI_low, CI_high, t,p))
}
demg = rbind(demotab2("age"),demg)
demg[2:(ncol(demg)-1)] = round(demg[2:(ncol(demg)-1)], 2)
demg$p = round(demg$p,4)
demg$CSIG = paste(demg$M_CSIG," (",demg$SD_CSIG,")",sep="")
demg$CLG = paste(demg$M_CLG," (",demg$SD_CLG,")",sep="")
demg$Difference = paste(demg$B," (",demg$SE,")",sep="")
demg$CI = paste(round(demg$CI_low,2), "; ",round(demg$CI_high,2),sep="")
demg = demg %>% dplyr::select(CLG,CSIG,Difference,CI,t,p)
write.csv(x = demg, paste(datapath,"PASAT_descriptive_continuous.csv",sep=""))
#
# finally, check median and mad for edss
median(unlist(na.omit(data%>%filter(CLG == 1))["edss"]))
mad(unlist(na.omit(data%>%filter(CLG == 1))["edss"]))
median(unlist(na.omit(data%>%filter(CLG == 0))["edss"]))
mad(unlist(na.omit(data%>%filter(CLG == 0))["edss"]))
# wilcox.test(edss~CLG,data)
# wilcox_effsize(data,edss~CLG,ci=T)
#
#
#
# Differences in frequencies ####
# CLG 1 = CLG
# sex: 1=female, 0=male
res=list()
for (i in 1:length(binaries)){
  # how the sausage is made, but obsolete for OR, CI, p:
  a = table(data.frame(data["CLG"],data[paste(binaries[i])]))[4] # both 1
  b = table(data.frame(data["CLG"],data[paste(binaries[i])]))[1] # both 0
  c = table(data.frame(data["CLG"],data[paste(binaries[i])]))[3] # only CLG 1
  d = table(data.frame(data["CLG"],data[paste(binaries[i])]))[2] # only var of interest 1
  # OR = (a * b) / (c * d) # resulting odds ratio
  # SE = sqrt(1/a + 1/b + 1/c + 1/d)
  # CI_low = exp(log(OR)-1.96*SE)
  # CI_high = exp(log(OR)+1.96*SE)
  # data.frame(OR,SE,CI_low, CI_high)
  res[[i]] = c(a,b,c,d,(fisher.test(table(data.frame(data["CLG"],data[paste(binaries[i])])))$estimate),
               unlist(fisher.test(table(data.frame(data["CLG"],data[paste(binaries[i])])))[2]),
               unlist(fisher.test(table(data.frame(data["CLG"],data[paste(binaries[i])])))[1]))
}
res = data.frame(list.rbind(res))
rownames(res) = binaries
res$OR=paste(round(res$odds.ratio,2)," [",round(res$conf.int1,2),",",round(res$conf.int2,2),"]",sep="")
res = res %>% dplyr::select(V1,V2,V3,V4,OR,p.value)
names(res) = c("CLG & Group membership = yes","CLG & Group membership = no", 
               "CSIG & Group membership = yes","CSIG & Group membership = no",
               "OR","p")
write.csv(x = res, paste(datapath,"PASAT_descriptive_binary.csv",sep=""))
#
#