# Logistic regression power/sensitivity analysis across models
# Max Korbmacher, 10.12.2024
#
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
# Predictors
predictors <- c("geno", "relapses_12mnths_before_baseline", "CH3L.1..mg.ml..mean", 
                "NfL..pg.ml.", "edss", "PASAT", "smoking_OFAMS", "BL_BMI", 
                "Omega3_suppl", "BAG_c", "baselineC", "baselineV", 
                "PF", "RF", "BP", "GH", "VT", "SF", "RE", "MH", "Vit_A_0", "Vit_D_0", 
                "Vit_E_0") # "Treatment_OFAMS", 
# Generate model formulas up to the 5th order to get an estimation of the probability of different events
max_order <- 8
all_combinations <- lapply(1:max_order, function(i) combn(predictors, i, simplify = FALSE))
all_formulas <- unlist(lapply(all_combinations, function(comb_list) {
  lapply(comb_list, function(comb) as.formula(paste("FLG ~ age +", paste(comb, collapse = " + "))))
}))
# Add optional terms (here, this is the covariate sex)
all_formulas <- c(all_formulas, lapply(all_formulas, function(f) update(f, . ~ . + sex)))
prob = function(form){
  a=glm(form,data=data,family = "binomial")
  p0=1/(1+exp(-1*(as.numeric(a$coefficients[1]))))
  p1=1/(1+exp(-1*(sum(a$coefficients))))
  return(c(p0,p1))
}
probdf=list()
for(i in 1:length(all_formulas)){
  probdf[[i]]=prob(all_formulas[[i]])
}
probdf = data.frame(list.rbind(probdf))
pwr = c()
for (i in 1:nrow(probdf)){
  pwr[i] = (wp.logistic(n = 88, p0 = (probdf$X1[i]), p1 = (probdf$X2[i]), alpha = 0.05,
                        power = NULL, alternative = c("greater"),
                        family = c("normal"), # we assume normal distribution of the predictor var@
                        parameter = NULL))$power
}
write.csv(x=pwr, file="power.csv")
