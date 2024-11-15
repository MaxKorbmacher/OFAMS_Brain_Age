# Multiverse to run all possible lms
print("Load packages and data.")
library(dplyr)
library(rlist)
data = read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
predictors = c("geno","relapses_12mnths_before_baseline",
               "CH3L.1..mg.ml..mean","NfL..pg.ml.",
               "edss_baseline" ,"smoking_OFAMS",
               "Mean_BMI_OFAMS","Treatment_OFAMS",
               "BAG_c","baselineC" ,"baselineV",
               "PF","RF","BP","GH","VT","SF","RE","MH")
print("Establishing all possible combinations of independently added variables to a list of formulas later used in LMs.")
# set up loop
vars = predictors
models = list()
for (i in 1:10){ # select all possible combinations up to the order of 10
  vc = combn(vars,i)
  for (j in 1:ncol(vc)){
    model <- as.formula(paste0("pasat_ARC~sex+Age_BL_OFAMS+", paste0(vc[,j], collapse = "+")))
    models <- c(models, model)
  }
  print(paste("Iteration ",i," done.",sep=""))
}
print("Adding blocks of variables to the formula list.")
# add vitamins as blocks (either present or not)
exp=expand.grid(models,c("","+Vit_A_0+Vit_D_0+Vit_E_0"))
full_vec = c()
for (i in 1:nrow(exp)){
  full_vec[i]=(paste(exp[i,1],exp[i,2],sep=""))
}
# add volume and control for TIV as block
exp=expand.grid(full_vec,c("","+TotalVol+TIV"))
full_vec = c()
for (i in 1:nrow(exp)){
  full_vec[i]=(paste(exp[i,1],exp[i,2],sep=""))
}
print("Running the regressions.")
# run regressions
coeflist=list()
for (i in 1:length(full_vec)){
  mods = lm(full_vec[i],data = data)
  coefs = data.frame(Beta = mods$coefficients,P=summary(mods)$coefficients[,4])
  #std_mods = lm(full_vec[i],data = data %>% mutate(across(where(is.numeric), scale)))
  #std_coefs=data.frame(Std.Beta = lm.beta(mods)$coefficients)
  #coeflist[[i]]=cbind(std_coefs,coefs)
  coeflist[[i]]=coefs
  coeflist[[i]]$Names = rownames(coeflist[[i]])
}
print("Extracting paramters for meta-analytic stats.")
# make a list of all predictors which we can then look at in terms of effect size / beta coefficient
allpreds=c(predictors,"sex","Age_BL_OFAMS", "Vit_A_0", "Vit_D_0", "Vit_E_0","TotalVol", "TIV")
res=temp=list()
for(j in 1:length(allpreds)){
  for(i in 1:length(coeflist)){
    temp[[i]]=coeflist[[i]] %>% filter(Names==allpreds[j])
  }
  res[[j]] = list.rbind(temp)
}
result = list.rbind(res)
print("Writing file.")
write.csv(result,"/Users/max/Documents/Local/MS/results/PASAT_effects.csv")
#
print("Done.")