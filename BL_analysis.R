# Baseline differences: functionally stable group (FSG) vs functional loss group (FLG)
# Max Korbmacher, 03.12.2024
#
# define data path
datapath = "/Users/max/Documents/Local/MS/results/"
# read data
data = read.csv(paste(datapath,"interrim_data.csv",sep=""))
# pkgs
library(dplyr)
#
# BL differemces between FLG and non-FLG patients ####
# Differences in continous scores ####

# age
FLG.imp = lm(age ~ FLG + sex, data)
summary(FLG.imp)

# Mean_BMI_OFAMS
FLG.imp = lm(Mean_BMI_OFAMS ~ FLG + age + sex, data)
summary(FLG.imp)

# EDSS
FLG.imp = lm(edss ~ FLG + age + sex,data)
summary(FLG.imp)
data %>% group_by(FLG) %>% summarize(M=mean(edss),SD=sd(edss))

# PASAT
FLG.imp = lm(PASAT ~ FLG + age + sex,data)
summary(FLG.imp)

# relapses_12mnths_before_baseline
FLG.imp = lm(relapses_12mnths_before_baseline ~ FLG + age + sex, data)
summary(FLG.imp)

# CH3L.1..mg.ml..mean
FLG.imp = lm(CH3L.1..mg.ml..mean ~ FLG + age + sex, data)
summary(FLG.imp)

# NfL..pg.ml.
FLG.imp = lm(NfL..pg.ml. ~ FLG + age + sex, data)
summary(FLG.imp)

# BAG_c
FLG.imp = lm(BAG_c ~ FLG + age + sex, data)
summary(FLG.imp)

# baselineC
FLG.imp = lm(baselineC ~ FLG + age + sex, data)
summary(FLG.imp)

# baselineV
FLG.imp = lm(baselineV ~ FLG + age + sex, data)
summary(FLG.imp)

# PF
FLG.imp = lm(PF ~ FLG + age + sex, data)
summary(FLG.imp)
max(na.omit(data$PF))
# RF
FLG.imp = lm(RF ~ FLG + age + sex, data)
summary(FLG.imp)

# BP
FLG.imp = lm(BP ~ FLG + age + sex, data)
summary(FLG.imp)

# GH
FLG.imp = lm(GH ~ FLG + age + sex, data)
summary(FLG.imp)

# VT
FLG.imp = lm(VT ~ FLG + age + sex, data)
summary(FLG.imp)

# SF
FLG.imp = lm(SF ~ FLG + age + sex, data)
summary(FLG.imp)

# RE
FLG.imp = lm(RE ~ FLG + age + sex, data)
summary(FLG.imp)

# MH                              
FLG.imp = lm(MH ~ FLG + age + sex, data)
summary(FLG.imp)



# Differences in frequencies ####
chisq.test((data$geno),data$FLG,simulate.p.value = TRUE, B = 1000000) # genotype diff
chisq.test((data$sex),data$FLG,simulate.p.value = TRUE, B = 1000000) # sex diff
chisq.test((data$smoking_OFAMS),data$FLG,simulate.p.value = TRUE, B = 1000000)
chisq.test((data$Treatment_OFAMS),data$FLG,simulate.p.value = TRUE, B = 1000000)
chisq.test((data$Omega3_suppl),data$FLG,simulate.p.value = TRUE, B = 1000000)


