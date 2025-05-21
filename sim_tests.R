# Analyse sim data
# includes checking VIF and some other things from fully saturated models
# Max Korbmacher, May 2025
#
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# libs
library(rlist)
library(dplyr)
library(regclass)
library(pROC)
library(stringr)
# 1. Simulated data -------
# dat
data = read.csv("/Users/max/Documents/Local/MS/results/synthetic_data.csv")
# mods
model = glm(FLG ~ .,
            data=data%>%select(-CLG,-RF,-RE),family = "binomial")
summary(model)
# vif and then exclusion
VIF(model)
# new model
model = glm(FLG ~ .,data=data%>%select(-CLG,-RF,-RE,-GH,-SF),family = "binomial")
VIF(model)
stats = data.frame(AUC = auc(model$data$FLG,predict(model,data)),Brier = DescTools::BrierScore(model),
           AIC = AIC(model),BIC = BIC(model))
# exclude SF-36
model1 = glm(FLG ~ .,data=data%>%select(-CLG,-RF,-RE,-GH,-SF,-PF,-BP,-VT,-MH),family = "binomial")
stats = rbind(stats,
      data.frame(AUC = auc(model1$data$FLG,predict(model1,data)),Brier = DescTools::BrierScore(model1),
           AIC = AIC(model1),BIC = BIC(model1)))
# removing MRI
model2 = glm(FLG ~ .,data=data%>%select(-CLG,-RF,-RE,-GH,-SF,-baselineC, -baselineV, -BAG_c),family = "binomial")

stats = rbind(stats,data.frame(AUC = auc(model2$data$FLG,predict(model2,data)),
                               Brier = DescTools::BrierScore(model2),
                               AIC = AIC(model2),BIC = BIC(model2)))
model3 = glm(FLG~.,data=data%>%select(FLG,edss, Vit_D_0, Vit_A_0,Current_DMT),family = "binomial")
data.frame(AUC = auc(model3$data$FLG,predict(model3,data)),Brier = DescTools::BrierScore(model3),
           AIC = AIC(model3),BIC = BIC(model3))
stats = rbind(stats,data.frame(AUC = auc(model3$data$FLG,predict(model3,data)),
                               Brier = DescTools::BrierScore(model3),
                               AIC = AIC(model3),BIC = BIC(model3)))

cbind( 	exp(coef(model3)), 	exp(summary(model3)$coefficients[,1] - 1.96*summary(model3)$coefficients[,2]), 	exp(summary(model3)$coefficients[,1] + 1.96*summary(model3)$coefficients[,2]) )


curious = data.frame(cbind( 	exp(coef(model)), 	exp(summary(model)$coefficients[,1] - 1.96*summary(model)$coefficients[,2]), 	exp(summary(model)$coefficients[,1] + 1.96*summary(model)$coefficients[,2]) ))
summary(model)$coefficients

l = data.frame(summary(model)$coefficients)
l = l[-1,]
l = l[-1,]
l$group = c("MRI","Demographic","Demographic",
            "Clinical","Clinical","Omics",
            "Clinical","Omics","Demographic",
            "Demographic","Intervention","MRI",
            "MRI", "SF-36","SF-36",
            "SF-36", "SF-36", "Omics",
            "Omics", "Omics","Clinical")

Names =row.names(l)
plot=ggforestplot::forestplot(
  df = l,
  name = Names,
  estimate = Estimate,
  se = Std..Error,
  xlab="Median odds ratio ± median absolute deviation",
  logodds = T,
  pvalue = Pr...z..,
  psignif = .05
)+
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
plot
#
# Note that the direction of the effect of EDSS makes no sense at all
# We exclude the variable
model = glm(FLG ~ .,data=data%>%select(-CLG,-RF,-RE,-GH,-SF,-edss),family = "binomial")
l = data.frame(summary(model)$coefficients)
l = l[-1,]
l$group = c("MRI","Demographic","Demographic",
            "Clinical","Clinical","Omics",
            "Clinical","Omics","Demographic",
            "Demographic", "Clinical","MRI",
            "MRI", "SF-36","SF-36",
            "SF-36", "SF-36", "Omics",
            "Omics", "Omics","Clinical")

Names =row.names(l)
plot=ggforestplot::forestplot(
  df = l,
  name = Names,
  estimate = Estimate,
  se = Std..Error,
  xlab="Median odds ratio ± median absolute deviation",
  logodds = T,
  pvalue = Pr...z..,
  psignif = .05
)+
  ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
plot






#exp(cbind(coef(model), confint(model)))  
#
#
#
model = glm(CLG ~ .,
            data=data%>%select(-FLG,-RF,-RE,-GH,-SF),family = "binomial")
summary(model)
cbind( 	exp(coef(model)), 	exp(summary(model)$coefficients[,1] - 1.96*summary(model)$coefficients[,2]), 	exp(summary(model)$coefficients[,1] + 1.96*summary(model)$coefficients[,2]) )

model = glm(CLG ~ .,data=data%>%select(-FLG,-RF,-RE,-GH,-SF),family = "binomial")
VIF(model)
stats = rbind(stats,data.frame(AUC = auc(model$data$CLG,predict(model,data)),Brier = DescTools::BrierScore(model),
                   AIC = AIC(model),BIC = BIC(model)))
# exclude SF-36
model1 = glm(CLG ~ .,data=data%>%select(-FLG,-RF,-RE,-GH,-SF,-PF,-BP,-VT,-MH),family = "binomial")
stats = rbind(stats,
              data.frame(AUC = auc(model1$data$CLG,predict(model1,data)),Brier = DescTools::BrierScore(model1),
                         AIC = AIC(model1),BIC = BIC(model1)))
# removing MRI
model2 = glm(CLG ~ .,data=data%>%select(-FLG,-RF,-RE,-GH,-SF,-baselineC, -baselineV, -BAG_c),family = "binomial")

stats = rbind(stats,data.frame(AUC = auc(model2$data$CLG,predict(model2,data)),
                               Brier = DescTools::BrierScore(model2),
                               AIC = AIC(model2),BIC = BIC(model2)))
model3 = glm(CLG~.,data=data%>%select(CLG,Current_DMT),family = "binomial")
data.frame(AUC = auc(model3$data$CLG,predict(model3,data)),Brier = DescTools::BrierScore(model3),
           AIC = AIC(model3),BIC = BIC(model3))
stats = rbind(stats,data.frame(AUC = auc(model3$data$CLG,predict(model3,data)),
                               Brier = DescTools::BrierScore(model3),
                               AIC = AIC(model3),BIC = BIC(model3)))

stats$Model = c("Saturated","SF-36 removed","MRI removed", "only significant")
stats$Prediction = c(replicate(nrow(stats)/2, "Disability Progression"), replicate(nrow(stats)/2, "Processing Speed Decline"))
stats[1:4] = round(stats[1:4],2)

write.csv(stats, "/Users/max/Documents/Local/MS/results/Simulation_Results.csv")




# 1.1 check model coefficients ----
#
# disability progression stable
curious$names = c(rownames(curious))
signi = c("edss","Current_DMT","age","Vit_A_0","Vit_D_0")
curious %>% filter(str_detect(as.character(names), str_c(as.character(signi), collapse="|")))
summary(model)
# processing speed stable
cbind( 	exp(coef(model)), 	exp(summary(model)$coefficients[,1] - 1.96*summary(model)$coefficients[,2]), 	exp(summary(model3)$coefficients[,1] + 1.96*summary(model)$coefficients[,2]) )




# Original data -------
data$FLG = factor(data$FLG)
data$FLG = relevel(data$FLG, ref = 2)
test = glm(FLG~.,data=data%>%dplyr::select(FLG,edss, Vit_D_0, Vit_A_0,Current_DMT),family = "binomial")
cbind( 	exp(coef(test)), 	exp(summary(test)$coefficients[,1] - 1.96*summary(test)$coefficients[,2]), 	exp(summary(test)$coefficients[,1] + 1.96*summary(test)$coefficients[,2]) )


