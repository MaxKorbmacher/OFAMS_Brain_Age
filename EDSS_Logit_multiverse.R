# Multiverse to run all reasonable Logistic Regression models for EDSS-based functional decline group classifications
#
# Max Korbmacher, 04.12.2024
#
# We use a GLM: model=glm(y_hat~x1+x2, family="binomial", data=data)
# And then we get both the coefficients and SEs
# These can be transformed to odds ratios with 95% CI: exp(coef(model)[i]±1.96*SE)
# Finally, also the 
#
# Note, the following assumptions informed the analytic choices:
# Collinearity
print("Load packages and data.")
library(dplyr)
library(rlist)
library(pscl)
library(AUC)
data = read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")
#
predictors = c("geno","relapses_12mnths_before_baseline",
               "CH3L.1..mg.ml..mean","NfL..pg.ml.",
               "edss","PASAT" ,"smoking_OFAMS","Mean_BMI_OFAMS","Treatment_OFAMS",
               "Omega3_suppl",
               "BAG_c","baselineC" ,"baselineV",
               "PF","RF","BP","GH","VT","SF","RE","MH",
               "Vit_A_0","Vit_D_0","Vit_E_0"
               )
print("Establishing all possible combinations of independently added variables to a list of formulas later used in LMs.")
# set up loop
vars = predictors
models = list()
for (i in 1:10){ # select all possible combinations up to the order of 10 (higher up, the model will be overparameterised & computation is not possible anymore)
  vc = combn(vars,i)
  for (j in 1:ncol(vc)){
    model <- as.formula(paste0("FLG~age+", paste0(vc[,j], collapse = "+")))
    models <- c(models, model)
  }
  print(paste("Iteration ",i," done.",sep=""))
}
# add sex or not
exp=expand.grid(models,c("","+sex"))
full_vec = c()
for (i in 1:nrow(exp)){
  full_vec[i]=(paste(exp[i,1],exp[i,2],sep=""))
}
# print("Adding blocks of variables to the formula list.")
# # add vitamins as blocks (either present or not)
# exp=expand.grid(models,c("","+Vit_A_0+Vit_D_0+Vit_E_0"))
# full_vec = c()
# for (i in 1:nrow(exp)){
#   full_vec[i]=(paste(exp[i,1],exp[i,2],sep=""))
# }
# # add volume and control for TIV as block
# exp=expand.grid(full_vec,c("","+TotalVol+TIV"))
# full_vec = c()
# for (i in 1:nrow(exp)){
#   full_vec[i]=(paste(exp[i,1],exp[i,2],sep=""))
# }
print("Running the regressions.")
# run regressions
coeflist=modinfo=list()
for (i in 1:length(full_vec)){
  mods = glm(full_vec[i],family="binomial",data = data)
  coefs = data.frame(Beta = mods$coefficients,SE = summary(mods)$coefficients[, 2], P=summary(mods)$coefficients[,4])
  #std_mods = lm(full_vec[i],data = data %>% mutate(across(where(is.numeric), scale)))
  #std_coefs=data.frame(Std.Beta = lm.beta(mods)$coefficients)
  #coeflist[[i]]=cbind(std_coefs,coefs)
  coeflist[[i]]=coefs
  coeflist[[i]]$Names = rownames(coeflist[[i]])
  modinfo[[i]] = c(pscl::pR2(mods)["McFadden"], full_vec[i]) # compute McFadden's pseudo R2
  #predicted = predict(mods, data, type="response") # predict
  #sensitivity(data$FLG, predicted) # can estimate sensitivity and specificity
  #specificity(data$FLG, predicted)
}
print("Overview of model evaluation criteria per model.")
modinfo = data.frame(list.rbind(modinfo))
print("Writing model metric table.")
write.csv(modinfo,"/Users/max/Documents/Local/MS/results/EDSS_model_info.csv")
#
print("Extracting parameters for meta-analytic stats.")
# make a list of all predictors which we can then look at in terms of effect size / beta coefficient
allpreds=c(predictors,"sex")
res=temp=list()
for(j in 1:length(allpreds)){
  for(i in 1:length(coeflist)){
    temp[[i]]=coeflist[[i]] %>% filter(Names==allpreds[j])
  }
  res[[j]] = list.rbind(temp)
}
result = list.rbind(res)
print("Writing file coefficients table.")
write.csv(result,"/Users/max/Documents/Local/MS/results/EDSS_stratification.csv")
#
print("Done.")
