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
#
data = data_unique(data = data, select = eid)

# 1. POWER ####
# power is calculated from another script (power.R or power_optimized.R for parallelisation)
# we use some descriptives on power and model fit
pwr = read.csv("/Users/max/Documents/Local/MS/results/model_info.csv")
# all models
hist(pwr$Power)
hist(pwr$McFaddenR2)
median(na.omit(pwr$Power))
mad(na.omit(pwr$Power))
median(na.omit(pwr$McFaddenR2))
mad(na.omit(pwr$McFaddenR2))
nrow(pwr) # number of all models ran
sum(is.na(pwr$Power)) # models with convergence errors disallowing power calculation
# high enough power models
pwr = pwr %>% filter(Power >= .8)
hist(pwr$Power)
hist(pwr$McFaddenR2)
median(na.omit(pwr$Power))
mad(na.omit(pwr$Power))
median(na.omit(pwr$McFaddenR2))
mad(na.omit(pwr$McFaddenR2))
nrow(pwr) # number of all models ran
# high enough power and excellent model fit
pwr = pwr %>% filter(McFaddenR2 >= .2)
hist(pwr$Power)
hist(pwr$McFaddenR2)
median(na.omit(pwr$Power))
mad(na.omit(pwr$Power))
median(na.omit(pwr$McFaddenR2))
mad(na.omit(pwr$McFaddenR2))
nrow(pwr) # number of all models ran
rm(pwr)
#
#
#
# 2. DEMOGRAPHICS: BL differences between FLG and non-FLG patients ####
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
contvar = c("BAG_c","baselineC","baselineV","edss","PASAT",
            "relapses_12mnths_before_baseline","age","BL_BMI",
            "CH3L.1..mg.ml..mean", "NfL..pg.ml.",
            "Vit_A_0", "Vit_D_0","Vit_E_0", "BP", "GH", "MH",
            "PF", "RE", "RF", "SF", "VT")
demg = list()
for (i in 1:length(contvar)){
  demg[[i]] = demotab(contvar[i])
}
demg = list.rbind(demg)
demg[2:(ncol(demg)-1)] = round(demg[2:(ncol(demg)-1)], 2)
demg$p = round(demg$p,4)
demg$FSG = paste(demg$M_FSG," (",demg$SD_FSG,")",sep="")
demg$FLG = paste(demg$M_FLG," (",demg$SD_FLG,")",sep="")
demg$Difference = paste(demg$B," (",demg$SE,")",sep="")
demg = demg %>% dplyr::select(FLG,FSG,Difference,CI_low,CI_high,t,p)
write.csv(x = demg, paste(datapath,"descriptive_continuous.csv",sep=""))
#
#
#
# Differences in frequencies ####
# FLG 1 = FLG
# sex: 1=female, 0=male
binaries = c("sex", "smoking_OFAMS","Omega3_suppl", "geno")
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