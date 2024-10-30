# train lifespan models from simple FS output
# then predict in several validation samples
# Max Korbmacher, 18 Oct 2024
#
# prep (load data and pkgs) ####
setwd("/cluster/projects/p33/users/maxk/MS/")
# load data
train=read.csv("training_set.csv") # HC training set
HC_validation=read.csv("HC_validation_set.csv") # HC validation set
long_validation=read.csv("long_validation_set.csv") # longitudinal validation set (ALZ & HC)
MS_validation=read.csv("MS_validation_set.csv")
# load pkgs
library(mgcv)
#install.packages("Metrics")
library(Metrics)
library(MASS)
#install.packages("relaimpo")
library(relaimpo)
library(lmtest)
library(lme4)
#library(sandwich)
#library(lm.beta)
library(effectsize)
#install.packages("rempsyc")
library(rempsyc)
library(dplyr)
library(ggeffects)
library(ggpubr)
library(marginaleffects)
library(MuMIn)
library(rlist)
library(lmerTest)
library(sjPlot)
#
# train ####
nms=train %>% dplyr::select(starts_with("lh_"), starts_with("rh_"))%>%names
#rhs=paste(nms,sep='',collapse = '+') # lm
#rhs=paste("age~",rhs,sep="",collapse = "+") # lm
rhs=paste('s(',nms,',k=4)',sep='',collapse = '+')
rhs=paste("age~",rhs,sep="",collapse = "+")
rhs=as.formula(rhs)
# model=gam(rhs, data = train)
# saveRDS(model, "results/model.rda")
model = readRDS("results/model.rda") # model can be loaded like this
#
# try also a linear model
rhs=paste(nms,sep='',collapse = '+') # lm
rhs=paste("age~",rhs,sep="",collapse = "+") # lm
rhs=as.formula(rhs)
model_lm=lm(rhs, data = train)
#
#
# explain models ####
# note that explanations of model parameters are quite limited in GAMs, due to the splines (difficult to compare these parameters)
params = parameters::model_parameters(model, robust = F) # unstandardized model coefficients
params_std = parameters::model_parameters(model, robust = T) # unstandardized model coefficients with robust SE (the same)
# params_std = parameters::model_parameters(model, standardize = "refit", robust = T, bootstrap = F, ci = 0.95) # standardized coefficients obtained by refitting (bootstrapping takes long, hence, we avoid it here).
std=data.frame(Parameter = nms, 
               df_robust = params_std$df[2:length(params_std$df)], p_robust = params_std$p[2:length(params_std$df)],
               F = params$`t / F`[2:length(params_std$df)], df_std = params$df[2:length(params_std$df)], p = params_std$p[2:length(params_std$df)])
write.csv(x = std, file = "results/explainability.csv",row.names = F)
# params = parameters::model_parameters(model, standardize = "basic", bootstrap = T, iterations = 1000, ci = 0.95) # This would be the bootstrapping alternative
# standardize_parameters(model, method = "posthoc", robust = T) # robust posthoc standardisation: scales parameters by median & MAD (instead of mean & SD)
#
#
#
#
# # other functions
# a=(lm.beta(model)) # simplest and equal to standardize_parameters() with method = "basic", only lm
# data.frame(a$standardized.coefficients)
# confint(model)
# coeftest(model,df=Inf)
# summary(model)
# coeftest(model, df = Inf, vcov = vcovHC) # in this large training sample = equal to standard settings in coeftest (no special covariance matrix)
# calc_rel = calc.relimp(model)
# boot_rel = boot.relimp(model, type = "lmg", nboot = 1000) # takes too long and eats too much memory
#
#
# predict in training sample
train$brainage = predict(model, train)
#
# retrieve intercept and slope from training predictions to obtain corrected brain age predictions
corlm = lm(brainage~age,data = train)
intercept = summary(corlm)$coefficients[1]
slope = summary(corlm)$coefficients[2]
correction_formula = "predicted_age_test + (age_test - (slope*predicted_age_test+intercept))"
corr_params = data.frame(intercept, slope, correction_formula)
write.csv(corr_params, "results/corr_params.csv", row.names = F)
correct = function(training_frame, predicted_age_train, predicted_age_test, age_test){
  corlm = lm(predicted_age_train~age,data = training_frame)
  intercept = summary(corlm)$coefficients[1]
  slope = summary(corlm)$coefficients[2]
  predicted_age_test + (age_test - (slope*predicted_age_test+intercept))
}
#
# evaluate
eval_metrics = function(model_type,age,pred.age){
  R=as.numeric(cor.test(age, pred.age)$estimate)
  R2=R^2
  MAE=mae(age, pred.age)
  RMSE=rmse(age, pred.age)
  return(data.frame(Model = model_type, R = R, R2 = R2, MAE=MAE, RMSE = RMSE))
}
train$brainage_corrected = correct(train,train$brainage,train$brainage,train$age)
train$brainage_corrected_lm = correct(train,predict(model_lm),predict(model_lm),train$age)

train_eval = rbind(eval_metrics("training: GAM", train$age, train$brainage),
                   eval_metrics("training: GAM corrected", train$age, train$brainage_corrected),
                   eval_metrics("training: LM", train$age, predict(model_lm)),
                   eval_metrics("training: LM corrected", train$age, train$brainage_corrected_lm)
      )
#
# validate model in different samples ####
# HC ####
# predict in test data with adapted eval formula
HC_validation$brainage = predict(model, HC_validation)
HC_validation$brainage_lm = predict(model_lm, HC_validation)
HC_validation$brainage_corrected = correct(train, train$brainage, HC_validation$brainage, HC_validation$age)
HC_validation$brainage_corrected_lm = correct(train, train$brainage, predict(model_lm, HC_validation), HC_validation$age)
HC_testeval = rbind(eval_metrics("HC validation: GAM", HC_validation$age, HC_validation$brainage),
      eval_metrics("HC validation: GAM corrected", HC_validation$age, HC_validation$brainage_corrected),
      eval_metrics("HC validation: LM", HC_validation$age, predict(model_lm, HC_validation)),
      eval_metrics("HC validation: LM corrected", HC_validation$age, HC_validation$brainage_corrected_lm)
)
train_res = rbind(train_eval, HC_testeval)
write.csv(train_res,"results/train_res.csv", row.names = F)
#hist(HC_validation$brainage_corrected_lm - HC_validation$age)
#hist(HC_validation$brainage_corrected - HC_validation$age)
#hist(HC_validation$brainage_corrected - HC_validation$age)
# Longitudinal HC, MCI, AD ####
# in longitudinal data
long_validation$subs = ifelse(grepl("sub-1",long_validation$eid),"sub-1",long_validation$eid) #label BBSC subs correctly
long_validation$subs = ifelse(grepl("sub-2",long_validation$subs),"sub-2",long_validation$subs)
long_validation$subs = ifelse(grepl("sub-3",long_validation$subs),"sub-3",long_validation$subs)
long_validation$brainage = predict(model, long_validation) # estimate brain ages
long_validation$brainage_lm=predict(model_lm, long_validation)
long_validation$brainage_corrected = correct(train, train$brainage, long_validation$brainage, long_validation$age)
long_validation$brainage_corrected_lm = correct(train, train$brainage, predict(model_lm, long_validation), long_validation$age)
ADNI=long_validation %>% filter(data=="ADNI")
ADNI_p1 = nice_scatter(
  data = ADNI,predictor = "age",  response = "brainage_corrected",
  ytitle = "Brain Age",  xtitle = "Age",
  has.points = FALSE,  has.jitter = F,  alpha = 1,
  has.confband = TRUE,  has.fullrange = FALSE,
  group = "diagnosis",  has.linetype = TRUE,
  has.shape = TRUE,  xmin = 75,
  xmax = 95,  xby = 5,  ymin = 75,  ymax = 95,
  yby = 5,  has.r = F,
  has.p = F,  #r.x = 5.5,  #r.y = 25,
  #colours = c("burlywood", "darkgoldenrod", "chocolate"),
  legend.title = "Diagnosis"
  #groups.labels = c("Weak", "Average", "Powerful")
)
ADNI_p1=ADNI_p1+theme(legend.position="bottom")+
  ggtitle("Within-subject variance naive ADNI ageing trajectories")+
  theme(plot.title = element_text(size = 14, face = "bold")) + 
  geom_abline(intercept = 0, slope = 1, size = 0.5) +
  theme(text = element_text(size=15))
ADNI_p2= nice_scatter( # each subject's slope
  data = ADNI,  predictor = "age",  response = "brainage_corrected",
  ytitle = "Brain Age",  xtitle = "Age", has.points = FALSE,
  has.jitter = F,  alpha = 1,
  has.confband = F,has.fullrange = FALSE,group = "subs",has.linetype = TRUE,
  has.shape = TRUE, xmin = 75,xmax = 95,xby = 5,ymin = 75,ymax = 95,yby = 5,
  has.legend = F
)
ADNI_p2=ADNI_p2+
  ggtitle("Diagnosis variance naive ADNI ageing trajectories")+
  theme(plot.title = element_text(size = 14, face = "bold"))+
  geom_abline(intercept = 0, slope = 1, size = 0.5)+
  theme(text = element_text(size=15))
#
# we could plot subs for each diagnosis, but that doesn't really make sense for 
# corrected predicted ages, as they appear very, very similar
# For uncorrected brain ages, these values are nearly unchanged / constant.
# Hence, we highly recommend to use corrected predicted age values!!!
# ADNI_p2.2= nice_scatter( # HC subjects' slopes
#   data = ADNI%>%filter(diagnosis=="AD"), predictor = "age",  response = "brainage_corrected",
#   ytitle = "Brain Age",  xtitle = "Age", has.points = FALSE,
#   has.jitter = T,  alpha = 1,
#   has.confband = F,has.fullrange = FALSE,group = "subs",has.linetype = TRUE,
#   has.shape = TRUE, xmin = 75,xmax = 100,xby = 5,ymin = 75,ymax = 100,yby = 5,
#   has.legend = F, has.r = T, has.p = T
# )
#
# arrange plots
ADNI_plot = egg::ggarrange(ADNI_p1,ADNI_p2, widths = c(1,1))
ggsave("results/ADNI_plot.pdf",ADNI_plot, width = 12, height = 7)
#ggsave("/tsd/p33/home/p33-maxk/ADNI_plot.pdf",ADNI_plot, width = 12, height = 7)
#
# check some demographics on the ADNI data
# test=(ADNI)%>%group_by(subs)%>%summarize(DIFF=max(age)-min(age))
# mean(test$DIFF)
# ADNI %>% filter(session=="ADNI Baseline") %>% summarize(MeanAge=mean(age),MinAge=min(age),MaxAge=max(age))
# 
# We show the Results in a table
ADNI$BAG=ADNI$brainage_corrected-ADNI$age
m0 = lmer(BAG~age+sex+diagnosis+(0+age|subs),data=ADNI)
m1 = lmer(BAG~age+sex+diagnosis+(1|diagnosis/subs),data=ADNI)
m2 = lmer(BAG~age+sex+diagnosis+(1|diagnosis),data=ADNI)
m3 = lmer(BAG~age+sex+diagnosis+(1|subs),data=ADNI)
r.squaredGLMM(m0) # the best model based on variances explained
r.squaredGLMM(m1)
r.squaredGLMM(m2)
r.squaredGLMM(m3)
anova(m0,m1,m2,m3) # m3 is however best/equal to m1 in terms of chisq test, aic, bic)
ggpredict(m0, terms=c("diagnosis","sex"), type = "re", interval = "prediction")
ggpredict(m1, terms=c("diagnosis","sex"), type = "re", interval = "prediction")
ggpredict(m2, terms=c("diagnosis","sex"), type = "re", interval = "prediction")
ggpredict(m3, terms=c("diagnosis","sex"), type = "re", interval = "prediction")

# we report both best performing models
Adj_BAG_ADNI = as.data.frame(rbind(
  ggpredict(m1, "diagnosis", type = "re", interval = "prediction"),
  ggpredict(m3, "diagnosis", type = "re", interval = "prediction"))
)
Adj_BAG_ADNI$group = c("lmer(BAG~age+sex+diagnosis+(1|diagnosis/subs)","lmer(BAG~age+sex+diagnosis+(1|diagnosis/subs)","lmer(BAG~age+sex+diagnosis+(1|diagnosis/subs)","lmer(BAG~age+sex+diagnosis+(1|subs)","lmer(BAG~age+sex+diagnosis+(1|subs)","lmer(BAG~age+sex+diagnosis+(1|subs)")
Adj_BAG_ADNI$age = "Age=85.45"
write.csv(Adj_BAG_ADNI,"results/Adj_BAG_ADNI.csv",row.names = F)
ggpredict(m2, terms=c("diagnosis", "age"), type = "re", interval = "prediction")

#
# The following are explorations and sanity checks of the brain age predictions.
# NOTE !!!! These predictions IGNORE within-subject variability.
# NOTE !!!! That means that data from the same subject which transitions between diagnoses is considered individually.
# In general, these reflect the marginal/predicted BAG from Adj_BAG_ADNI.
# Namely, brain age is (generally) more accurately predicted in HC than in MCI and AD
#
#
TrainMetr=rbind(
#Linear Models
  cbind(long_validation %>% group_by(diagnosis) %>% 
          summarize(model="LM_raw",R=cor(age, brainage_lm),
                    MAE=mae(age, brainage_lm), RMSE=rmse(age, brainage_lm)),
        long_validation %>% group_by(diagnosis) %>% 
          group_map(~ broom::tidy(lm(brainage_lm ~ age+sex, data = .x))[2,]) %>% list.rbind),
  cbind(long_validation %>% group_by(diagnosis) %>% 
          summarize(model="LM_cor",R=cor(age, brainage_corrected_lm),
                    MAE=mae(age, brainage_corrected_lm), RMSE=rmse(age, brainage_corrected_lm)),
        long_validation %>% group_by(diagnosis) %>% 
          group_map(~ broom::tidy(lm(brainage_corrected_lm ~ age+sex, data = .x))[2,]) %>% list.rbind),
# GAMs
  cbind(long_validation %>% group_by(diagnosis) %>% 
          summarize(model="GAM_raw",R=cor(age, brainage),
                    MAE=mae(age, brainage), RMSE=rmse(age, brainage)),
        long_validation %>% group_by(diagnosis) %>% 
          group_map(~ broom::tidy(lm(brainage ~ age+sex, data = .x))[2,]) %>% list.rbind),
  cbind(long_validation %>% group_by(diagnosis) %>% 
          summarize(model="GAM_cor",R=cor(age, brainage_corrected),
                    MAE=mae(age, brainage_corrected), RMSE=rmse(age, brainage_corrected)),
        long_validation %>% group_by(diagnosis) %>% 
          group_map(~ broom::tidy(lm(brainage_corrected ~ age+sex, data = .x))[2,]) %>% list.rbind)
)
write.csv(TrainMetr,"results/ADNI_group_metrics.csv")
#
#
# Another way to sanity check the data, potentially useful later
# rbind(long_validation %>% filter(diagnosis=="HC") %>% summarize(diagnosis="HC",model="LM_raw",R=cor(age, brainage_lm),MAE=mae(age, brainage_lm), RMSE=rmse(age, brainage_lm)), #print("Linear model BAG.")
#       long_validation %>% filter(diagnosis=="MCI") %>% summarize(diagnosis="MCI",model="LM_raw",R=cor(age, brainage_lm),MAE=mae(age, brainage_lm), RMSE=rmse(age, brainage_lm)),
#       long_validation %>% filter(diagnosis=="AD") %>% summarize(diagnosis="AD",model="LM_raw",R=cor(age, brainage_lm),MAE=mae(age, brainage_lm), RMSE=rmse(age, brainage_lm)),
#       long_validation %>% filter(diagnosis=="HC") %>% summarize(diagnosis="HC",model="LM_cor",R=cor(age, brainage_corrected_lm),MAE=mae(age, brainage_corrected_lm), RMSE=rmse(age, brainage_corrected_lm)), #print("Corrected linear model BAG.")
#       long_validation %>% filter(diagnosis=="MCI") %>% summarize(diagnosis="MCI",model="LM_cor",R=cor(age, brainage_corrected_lm),MAE=mae(age, brainage_corrected_lm), RMSE=rmse(age, brainage_corrected_lm)),
#       long_validation %>% filter(diagnosis=="AD") %>% summarize(diagnosis="AD",model="LM_cor",R=cor(age, brainage_corrected_lm),MAE=mae(age, brainage_corrected_lm), RMSE=rmse(age, brainage_corrected_lm)),
#       long_validation %>% filter(diagnosis=="HC") %>% summarize(diagnosis="HC",model="GAM_raw",R=cor(age, brainage),MAE=mae(age, brainage), RMSE=rmse(age, brainage)), #
#       long_validation %>% filter(diagnosis=="MCI") %>% summarize(diagnosis="MCI",model="GAM_raw",R=cor(age, brainage),MAE=mae(age, brainage), RMSE=rmse(age, brainage)),
#       long_validation %>% filter(diagnosis=="AD") %>% summarize(diagnosis="AD",model="GAM_raw",R=cor(age, brainage),MAE=mae(age, brainage), RMSE=rmse(age, brainage)), #print("Corrected GAM BAG.")
#       long_validation %>% filter(diagnosis=="HC") %>% summarize(diagnosis="HC",model="GAM_cor",R=cor(age, brainage_corrected),MAE=mae(age, brainage_corrected), RMSE=rmse(age, brainage_corrected)),
#       long_validation %>% filter(diagnosis=="MCI") %>% summarize(diagnosis="MCI",model="GAM_cor",R=cor(age, brainage_corrected),MAE=mae(age, brainage_corrected), RMSE=rmse(age, brainage_corrected)),
#       long_validation %>% filter(diagnosis=="AD") %>% summarize(diagnosis="AD",model="GAM_cor",R=cor(age, brainage_corrected),MAE=mae(age, brainage_corrected), RMSE=rmse(age, brainage_corrected)))
# single subject sanity checks
#test = long_validation %>% group_by(subs,diagnosis) %>% summarize(R=cor(age, brainage), MAE=mae(age, brainage_corrected_lm), RMSE=rmse(age, brainage))
#test = long_validation %>% group_by(subs,diagnosis) %>% summarize(R=cor(age, brainage_corrected), MAE=mae(age, brainage_corrected), RMSE=rmse(age, brainage_corrected))
# long_validation %>% filter(subs=="021_S_4659") %>% dplyr::select(age, brainage,brainage_corrected, brainage_lm, brainage_corrected_lm)
# long_validation %>% filter(subs=="011_S_0861") %>% dplyr::select(age, brainage,brainage_corrected, brainage_lm, brainage_corrected_lm)
#
#
# BBSC
BBSC=long_validation%>%filter(data=="BBSC")
#BBSC[BBSC$session %in%  (c("sub-1_1","sub-2_1","sub-3_1")),] %>% select(age)
#nrow(BBSC)
BBSC_p1 = nice_scatter(
  data = BBSC,predictor = "age",  response = "brainage_corrected",
  ytitle = "Brain Age",  xtitle = "Age",
  has.points = FALSE,  has.jitter = T,  alpha = 1,
  has.confband = TRUE,  has.fullrange = FALSE,
  group = "subs",  has.linetype = TRUE, has.shape = TRUE,  xmin = 25,
  xmax = 45,  xby = 5,  ymin = 25,  ymax = 45, yby = 5,  has.r = F,
  has.p = F,  #r.x = 5.5,  #r.y = 25,
  colours = c("burlywood", "darkgoldenrod", "chocolate"),
  # legend.title = "Diagnosis"
  #groups.labels = c("Weak", "Average", "Powerful")
)
BBSC_p1 = BBSC_p1+geom_abline(intercept = 0, slope = 1, size = 0.5)+
  ggtitle("Within-subject variance naive BBSC ageing trajectories")+
  theme(plot.title = element_text(size = 14, face = "bold"))+
  theme(text = element_text(size=15))
m4=lmer(brainage_corrected~age+(1|subs),data=BBSC)
summary(m4)
BBSC%>%group_by(subs)%>%summarize(R=cor(brainage_corrected,age), R_raw=cor(brainage,age))

# Model testing (MS vs HC) ####
MS_validation$brainage=predict(model, MS_validation)
MS_validation$brainage_lm=predict(model_lm, MS_validation)
MS_validation$brainage_corrected = correct(train, train$brainage, MS_validation$brainage, MS_validation$age)
MS_validation$brainage_corrected_lm = correct(train, train$brainage, predict(model_lm, MS_validation), MS_validation$age)
MS_validation$BAG=MS_validation$brainage_corrected-MS_validation$age
MS=MS_validation%>%filter(diagnosis=="MS")
MS_control=MS_validation%>%filter(diagnosis=="HC")
MS_validation_stats=rbind(eval_metrics("MS:GAM", MS$age, MS$brainage),
      eval_metrics("MS:GAM corrected", MS$age, MS$brainage_corrected),
      eval_metrics("MS:LM", MS$age, MS$brainage_lm),
      eval_metrics("MS:LM corrected", MS$age, MS$brainage_corrected_lm),
      
      eval_metrics("HC:GAM", MS_control$age, MS_control$brainage),
      eval_metrics("HC:GAM corrected", MS_control$age, MS_control$brainage_corrected),
      eval_metrics("HC:LM", MS_control$age, MS_control$brainage_lm),
      eval_metrics("HC:LM corrected", MS_control$age, MS_control$brainage_corrected_lm)
)
write.csv(MS_validation_stats,"results/MS_validation_stats.csv",row.names = F)
# raw differences
cohens_d(MS_control$BAG,MS$BAG)
test=t.test(MS_control$BAG,MS$BAG)

cohens_d((MS_control%>%filter(sex=="F"))$BAG,(MS%>%filter(sex=="F"))$BAG)
cohens_d((MS_control%>%filter(sex=="M"))$BAG,(MS%>%filter(sex=="M"))$BAG)

# check also for scanner (RE), sex, age as covariates
m5 = lmer(BAG~diagnosis+age+sex+(1|scanner),data=MS_validation)
m6 = lm(BAG~diagnosis+age+sex+scanner,data=MS_validation)
summary(m5)
r.squaredGLMM(m5)
summary(m6)
anova(m5,m6) # the linear model performs better and will be used
#ggpredict(m6, terms=c("diagnosis","scanner"), type = "fe")
#ggpredict(m6, terms=c("diagnosis","sex"), type = "fe")
m7 = lm(BAG~diagnosis*sex+age+scanner,data=MS_validation)
summary(m7)
avg_comparisons(m5,variables = c("diagnosis","sex","age"))
avg_comparisons(m6,variables = c("diagnosis","sex","age"))
p1=plot_comparisons(m7, variables = "diagnosis", by="sex", comparison="difference")
p1=p1+theme_classic()+ylab("BAG Difference in years, 95% CI")+xlab("Sex")+
  theme(text = element_text(size=15)) + ggtitle("Sex-specific BAG differences: MS vs HC")+
  theme(plot.title = element_text(size = 14, face = "bold"))
#test=predictions(m5)
#print(test, style="data.frame")
ggpredict(m6, terms=c("diagnosis","sex"), type = "fe")
#
#
# put all plots together
val_plot=egg::ggarrange(p1,BBSC_p1, ADNI_p1,ADNI_p2, ncol=2,nrow=2, widths = c(1,1),heights = c(1,1))
ggsave(filename = "results/val_plot.pdf",val_plot, width = 12, height = 10)
#ggsave("/tsd/p33/home/p33-maxk/val_plot.pdf",val_plot, width = 12, height = 10)
#
#
#
print("The GAM outperforms the (simpler) linear model using raw brain age values.")
print("The models are similar after training sample age bias correction.")
print("Still, the GAM-based corrected brain age gaps differentiate AD from MCI (and HC) better than the linear models.")
#
#
# BAG longitudinal checks ####
# for comparability with the target data, we are interested in expectable trajectories
# here, we have BBSC and ADNI data available
sjPlot::plot_model(m3, type="pred", terms=c("age","diagnosis","sex"), pred.type="fe", ci.lvl=.95)
p1=ggplot(data= ADNI%>%filter(diagnosis=="HC"), aes(x = age, y = BAG)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  stat_smooth(method = "lm",aes(group = 1),color="blue",fill = "blue")+ #
  ylab("BAG HC") + xlab("Age") +
  #stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
p2=ggplot(data= ADNI%>%filter(diagnosis=="MCI"), aes(x = age, y = BAG)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  stat_smooth(method = "lm",aes(group = 1),color="blue",fill = "blue")+ #
  ylab("BAG MCI") + xlab("Age") +
  #stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
p3=ggplot(data= ADNI%>%filter(diagnosis=="AD"), aes(x = age, y = BAG)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  stat_smooth(method = "lm",aes(group = 1),color="blue",fill = "blue")+ #
  ylab("BAG AD") + xlab("Age") +
  #stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
ADNI_plot2=ggarrange(p1,p2,p3)
ggsave("results/ADNI_plot2.pdf",ADNI_plot2, width = 6, height = 6)
ggsave("/tsd/p33/home/p33-maxk/ADNI_plot2.pdf",ADNI_plot2, width = 6, height = 6)

long_validation$BAG=long_validation$brainage_corrected-long_validation$age

p1=ggplot(data= long_validation%>%filter(data=="BBSC")%>%filter(subs=="sub-1"), aes(x = age, y = BAG)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = subs), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  stat_smooth(method = "lm",aes(group = 1),color="blue",fill = "blue")+ #
  ylab("BAG BBSC1") + xlab("Age") +
  #stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
p2=ggplot(data= long_validation%>%filter(data=="BBSC")%>%filter(subs=="sub-2"), aes(x = age, y = BAG)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = subs), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  stat_smooth(method = "lm",aes(group = 1),color="blue",fill = "blue")+ #
  ylab("BAG BBSC2") + xlab("Age") +
  #stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
p3=ggplot(data= long_validation%>%filter(data=="BBSC")%>%filter(subs=="sub-3"), aes(x = age, y = BAG)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = subs), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  stat_smooth(method = "lm",aes(group = 1),color="blue",fill = "blue")+ #
  ylab("BAG BBSC3") + xlab("Age") +
  #stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
BBSC_p2=ggarrange(p1,p2,p3)
ggsave("results/BBSC_p2.pdf",BBSC_p2, width = 6, height = 6)
ggsave("/tsd/p33/home/p33-maxk/BBSC_p2.pdf",BBSC_p2, width = 6, height = 6)
hist(train$age,breaks = 100)
nrow(train)
