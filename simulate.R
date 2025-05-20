# Synth MS
# Max Korbmacher
# May 2025
#
#install.packages("synthpop")
library(synthpop)
library(rlist)
data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
data = data %>% select("CLG","Current_DMT","geno", "relapses_12mnths_before_baseline", 
                                "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI",
                                "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
                                "Vit_A_0","Vit_D_0", "Vit_E_0",
                                "PF", "BP", "GH", "VT", "SF","MH","sex","age")
mysyn = syn(data,m=10)
summary(mysyn)
data = list.rbind(mysyn$syn)

model = glm(CLG ~ age + Current_DMT + sex + geno+ relapses_12mnths_before_baseline+ 
                edss+ PASAT+ smoking_OFAMS+
                Vit_A_0+Vit_D_0+ Vit_E_0,
                data=data,family = "binomial")
summary(model)

#
#
#
#
#
#


data <- read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
data = data %>% select("FLG","Current_DMT","geno", "relapses_12mnths_before_baseline", 
                       "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI",
                       "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
                       "Vit_A_0","Vit_D_0", "Vit_E_0",
                       "PF", "BP", "GH", "VT", "SF","MH","sex","age")
mysyn = syn(data,m=10)
summary(mysyn)
data = list.rbind(mysyn$syn)

model = glm(FLG ~ age + Current_DMT + sex + geno+ relapses_12mnths_before_baseline+ 
              edss+ PASAT+ smoking_OFAMS+
              Vit_A_0+Vit_D_0+ Vit_E_0,
            data=data,family = "binomial")
summary(model)






install.packages("finalfit")
library("finalfit")
library(dplyr)
explanatory = c("Current_DMT","geno", "relapses_12mnths_before_baseline", 
  "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI",
  "Omega3_suppl", "BAG_c", "baselineC", "baselineV",
  "Vit_A_0","Vit_D_0", "Vit_E_0",
  "PF", "BP", "GH", "VT", "SF","MH","sex","age")
dependent="CLG"

data %>%
  finalfit_newdata(explanatory = explanatory, newdata = list(
    c("<40 years",  "Submucosa", "No"),
    c("<40 years", "Submucosa", "Yes")
    )) -> newdata

data %>% glmmulti(dependent, explanatory) %>% 
  boot_predict(newdata, 
               estimate_name = "Predicted probability of death",
               R=100, boot_compare = FALSE,
               digits = c(2,3))
