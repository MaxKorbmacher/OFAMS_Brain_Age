# Synth MS
# Max Korbmacher
# May 2025
#
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
#
# packages
# install.packages("synthpop")
# install.packages("regclass")
library(synthpop)
library(rlist)
library(dplyr)
library(regclass)
#
# data
data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
names(data)
data = data%>%select(-X,-CH3L.1..mg.ml..mean,-SYMPTOM_DURATION,-DISEASE_DURATION, -age_gap,-eid, -Omega3_suppl)
# data = data %>% select("CLG","Current_DMT","geno", "relapses_12mnths_before_baseline", 
#                                 "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI",
#                                 "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
#                                 "Vit_A_0","Vit_D_0", "Vit_E_0",
#                                 "PF", "BP", "GH", "VT", "SF","MH","sex","age")
# factor specification
data$FLG = factor(data$FLG)
data$CLG = factor(data$CLG)
data$sex = factor(data$sex)
data$smoking_OFAMS = factor(data$smoking_OFAMS)
data$Treatment_OFAMS = factor(data$Treatment_OFAMS)
data$Current_DMT = factor(data$Current_DMT)
# creating synthetic data
mysyn = syn(data,m=100,k=10000,seed = 123)
summary(mysyn)
data = list.rbind(mysyn$syn)
write.csv(data, "/Users/max/Documents/Local/MS/results/synthetic_data.csv")
