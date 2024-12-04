# OFAMS brain age analysis
# 23 October 2024
# Max Korbmacher
#
# Note: this is my messy analysis file to clean/wrangle and get a first impression of the data and plot some trends, etc.
# The most important is the produced interrim_data.csv which provides the input for the multiverse analyses.
#
# set savepath
savepath = "/Users/max/Documents/Local/MS/results/"
# 0.1 Prep ####
# load pkgs
# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(haven,car,dplyr,lme4,lmerTest,ggeffects,sjPlot,
               data.table,factoextra,simputation,ggpubr,psych,reshape2,
               gplots,rstatix,ggseg,marginaleffects,MuMIn,tidyr,#glmnet,ipflasso,parameters,
               olsrr,parallel,gghighlight,cowplot,zoo)
#
#
# load data
data = read.csv("/Users/max/Documents/Local/MS/GAM_preds.csv") # contains brain age predictions
edss = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/edss.sav") # contains Expanded Disability Status Scale (EDSS) scores
prom = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/sf36.sas7bdat") # patient reported outcome measures (PROM)
fati = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/fatigue.sas7bdat") # fatigue scores
demo = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/Demographics.sav") # Demographics
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
relapse = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/relapsereport.sas7bdat") # relapse info
geno = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/hla_drb1_15_Juni_2011_alle pas.sav") # genotype 
icv = read.csv("/Users/max/Documents/Local/MS/data/icv.csv")
lesion_count = read.csv("/Users/max/Documents/Local/MS/data/lesion_count.csv")
lesion_vol = read.csv("/Users/max/Documents/Local/MS/data/lesion_vol.csv")
blood = read.csv("/Users/max/Documents/Local/MS/data/Copy of Paired_serum-MRI.csv")
#
# open questions about data:
# should any inflammation markers be included? E.g. inflamm88.sav (or some of the files I don't understand?)
# files I don't understand: fa_bl.sas7bdat, hla_drb1_15_Juni_2011_alle pas.sav 
# finally, here are the baseline measures:
# demo = read_spss("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/BASELINE.sav") # baseline info
#
#
# simple overview
# data %>% group_by(session) %>% summarize(Mage = mean(age),
#                                          M_corrected_BAG = mean(corrected_brainage-age),
#                                          M_BAG = mean(brainage-age), N = length(brainage))
# cor(data$corrected_brainage-data$age, data$session)
# cor(data$corrected_brainage,data$age)
# cor(data$brainage,data$age)
#
# copy original data frame for later analysis of brain structure 
anat = data
# caluculate BAG
data = data %>% select(eid, session, age, brainage, corrected_brainage)
data$BAG_c = data$corrected_brainage-data$age
data$BAG = data$brainage-data$age
#
#
######## ADD COVARIATES!
# start with edss
#edss %>% select(Patnr,contains("EDSS")) %>% melt(id.vars="Patnr")
edss_long = melt(edss, id.vars="Patnr")
names(edss_long) = c("eid","session","edss")
print("Note that we only have EDSS scores for baseline and certain follow-up times.")
print("We can note an increase in EDSS scores over time. (Yet, with increasing variablility.)")
edss_long %>% group_by(session) %>% summarize(M=mean(na.omit(edss)), SD=sd(na.omit(edss))) # show mean+sd for each session
levels(edss_long$session) = c(0,6,12,18,24) # standardize session names
data$session = factor(data$session) # treat session as the factor it is also in the brain age data frame
data = full_join(data,edss_long, by=c("eid","session")) # join edss and brain age data frames.
#
# demographics: demo$Gender # 0 = Male
# genotype: geno$HLA_1501_1 # HLADRB1501 positive = 1
# relapse: relapse$relapseno # nb of relapses > var of interest
demo = right_join(demo, geno, by = "Patnr")
relapse$Patnr = relapse$patno
demo = right_join(demo, relapse, by = "Patnr")
print("This completes the cross-sectional / baseline data in the demo frame.")
# patient recorded outcome measures (prom)
prom %>% select(patno,contains("spm")) %>% is.na() %>% colSums() # check number of NAs (none)
#print("There are some (very few) NAs. We impute using sequential hot deck imputation.")
#prom=prom %>% select(patno,VISIT,contains("spm"))
#prom=impute_shd(prom, .~1, backend="VIM")
visits = levels(factor(prom$VISIT))
# old code for pca (if desired)
# prom_pca = prom_pca_plot = list()
# for (i in 1:length(visits)){
#   prom_pca[[i]] = prom %>% filter(VISIT == visits[i]) %>% select(contains("spm")) %>% prcomp
#   prom_pca_plot[[i]] = fviz_eig(prom_pca[[i]], main = paste("PROM PCA ",visits[i],sep=""), addlabels=TRUE, hjust = -0.3,linecolor ="red") + theme_minimal() + ylim(0,100)
# }
# ggarrange(plotlist = prom_pca_plot) # visualise the variance explained by each factor
# prom$PROM_PC = c(prom_pca[[1]]$x[,1],prom_pca[[2]]$x[,1],prom_pca[[3]]$x[,1],prom_pca[[4]]$x[,1])# First factor explains > 37%. We can use that one.
prom = prom %>% rename(eid = patno, session = VISIT)
prom$session=factor(prom$session)
levels(prom$session)=c(0,12,24,6)
data=right_join(data,prom,by=c("eid","session"))
#
# Add sex and genotype to demo data.
geno = geno %>% rename(eid=Patnr,geno=HLA_1501_1)
sex = demo %>% select(Patnr,Gender) %>% rename(eid=Patnr,sex=Gender)
data = right_join(data,geno, by="eid")
data = right_join(data,sex, by="eid")
#
#
#
#
# 0.2 Functions ####
simple_range_extracter <- function(p, scale) {
  d <- ggplot_build(p)
  d$plot$scales$get_scales(scale$aesthetics)$range$range
}
get_shared_scale <- function(..., scale) {
  plots <- list(...)
  ranges <- purrr::map(plots, ~simple_range_extracter(., scale))
  single_range <- range(unlist(ranges))
  scale$limits <- single_range
  scale
}
# Main function
set_scale_union <- function(..., scale) {
  exprs <- rlang::enexprs(...)
  scale <- get_shared_scale(..., scale = scale)
  var_nms <- purrr::map_chr(exprs, rlang::as_name)
  edit_plots_in_place(var_nms, env = parent.frame(),
                      scale = scale)
  # Invisibly return the scale, in case you need it later
  invisible(scale)
}

# Sub-function
edit_plots_in_place <- function(names, env, scale) {
  vars <- rlang::env_has(env = env, nms = names)
  if (!all(vars))
    stop("Environment does not have variables for ",
         paste(names(vars[!vars]), collapse=", "))
  purrr:::walk(names, function(nm) {
    og_plot <- rlang::env_get(env, nm = nm)
    message("Changing plot `", nm, "`")
    # Muffles messages about already having scales
    withCallingHandlers(
      assign(x = nm, envir = env,
             value = og_plot + scale),
      message = function(err) {
        if (grepl("already present", err$message))
          invokeRestart("muffleMessage")
      })
  })
}
# 1.1 Descriptives ####

# correlation matrix of proms
x = cor(data%>%select(contains("spm")),use = "complete.obs")
col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = x, col = col, symm = TRUE)
print("The heatmap underscores the findings by Hobart et al. (2001) https://doi.org/10.1136/jnnp.71.3.363 that the subscale sum scores are not that straight forward.")
print("We will anyways go ahead with summing into 8 categories.")
# add proms sum scores (for interpretability divided by nb of item per category)
# note that the scale levels differ and absolute sum scores have to be used.
data$PF = data%>%select(contains("spm3")) %>% rowSums() # Physical functioning (PF)
data$RF = data%>%select(contains("spm4")) %>% rowSums() # Role-physical (RF)
data$BP = data%>%select(spm7,spm8) %>% rowSums() # Bodily pain (BP)
data$GH = data%>%select(spm1,contains("spm11")) %>% rowSums() # General health (GH)
data$VT = data%>%select(spm9a,spm9b,spm9c,spm9d) %>% rowSums() # Vitality (VT)
data$SF = data%>%select(spm6,spm10) %>% rowSums() # Social functioning (SF)
data$RE = data%>%select(contains("spm5")) %>% rowSums() # Role-emotional (RE)
data$MH = data%>%select(spm9b, spm9c, spm9d, spm9f, spm3h) %>% rowSums() # Mental health (MH)
#
# correlation matrices of variables of interest
visits = c("0","6","12","24") # at 120 months follow-up, there are no clinical measures taken.
df = unique(data)
for (i in 1:length(visits)){
  x = cor(df%>% filter(session == visits[i])%>%select(age,BAG,BAG_c,edss,PF,RF,BP,GH,VT,SF,RE,MH),use = "complete.obs")
  pdf(paste(savepath,"heatmap_",visits[i],"_months",".pdf",sep=""), width=8, height=8)
  heatmap.2(Rowv = F, x = x, col = col, symm = T, cellnote=round(x,2),notecol="black", dendrogram = "none",key=T,density.info="density",trace = "none",notecex = 1.25)
  dev.off()
}
# since there are only very few brain age values at 12 months, we make a cor matrix without it
x = cor(df%>% filter(session == visits[i])%>%select(edss,PF,RF,BP,GH,VT,SF,RE,MH),use = "complete.obs")
pdf(paste(savepath,"heatmap_",visits[i],"_months",".pdf",sep=""), width=8, height=8)
heatmap.2(Rowv = F, x = x, col = col, symm = T, cellnote=round(x,2),notecol="black", dendrogram = "none",key=T,density.info="density",trace = "none",notecex = 1.25)
dev.off()
#
print("There are multiple duplicates in the dataset")
print(paste("To be exact: ",nrow(data)-nrow(unique(data))," rows.",sep=""))
print("Another problem are the spread of participants in different session bins:")
data %>% group_by(session) %>% summarize(N=length(na.omit((BAG))))
print("We combat these problems in the representation of initial descriptives by excluding duplicates and low-N sessions.")
df2 = filter(df,session == 0 | session == 6 | session == 24 |session == 120)
print("All time points.")
ggplot(df, aes(session, BAG_c, group = eid, color = ""))+ #geom_point() + 
  geom_line() + scale_colour_manual(values = c("brown")) +
  theme_classic() + theme(legend.position = "none")
 print("Only time points with sufficient data.")
ggplot(df2, aes(session, BAG_c, group = eid, color = ""))+ #geom_point() + 
   geom_line() + scale_colour_manual(values = c("brown")) +
   theme_classic() + theme(legend.position = "none")
print("All time points.")
ggplot(data= df, aes(x = session, y = BAG_c)) + 
  geom_point(size = 2, alpha= 1, color = "gray85") + #,  aes(color = eid) colour points by group
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="blue",fill = "red")+ #method = "lm"
  ylab("Brain Age Gap") + xlab("Session") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
print("Only time points with sufficient data.")
BAG_descriptive = ggplot(data= df2, aes(x = session, y = BAG_c)) + 
  #scale_shape_manual(values=1:nlevels(factor(df2$eid))) +
  geom_point(size = 1.5, alpha= 1, color = "gray85") + #,  aes(color = eid) colour points by group # shape by eid aes(shape=factor(eid)),
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Brain Age Gap (years)") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
ggsave(paste(savepath,"BAG_descriptive.pdf",sep=""),BAG_descriptive)
df%>%group_by(session)%>%summarize(M=mean(na.omit(edss)))
df3 = filter(df,session == 0 | session == 6 | session == 12 |session == 18 |session == 24)
level_order = c(0,6,12,18,24)
# EDSS_descriptive = ggplot(data= df3, aes(x = factor(session,level=level_order), y = edss)) + 
#   geom_point(size = 1.5, alpha= 1, color = "gray85") +
#   geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
#   stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
#   ylab("EDSS") + xlab("Session (months)") +
#   stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
#   theme_bw()
# #setwd(paste(savepath,sep=""))
df4 = filter(df,session == 0 | session == 6 | session == 12 |session == 24)
PF_descriptive = ggplot(data= df4, aes(x = session, y = PF)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Physical functioning") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
RF_descriptive = ggplot(data= df4, aes(x = session, y = RF)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Role-physical") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
BP_descriptive = ggplot(data= df4, aes(x = session, y = BP)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Bodily pain") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
GH_descriptive = ggplot(data= df4, aes(x = session, y = GH)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("General health") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
Physical_Health = ggarrange(PF_descriptive,RF_descriptive,BP_descriptive,GH_descriptive)
ggsave(paste(savepath,"Physical_Health.pdf",sep=""),Physical_Health)
VT_descriptive = ggplot(data= df4, aes(x = session, y = VT)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("General health") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
SF_descriptive = ggplot(data= df4, aes(x = session, y = SF)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Social functioning") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
RE_descriptive = ggplot(data= df4, aes(x = session, y = RE)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Role-emotional") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
MH_descriptive = ggplot(data= df4, aes(x = session, y = MH)) + 
  geom_point(size = 1.5, alpha= 1, color = "gray85") +
  geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
  stat_smooth(aes(group = 1),color="red",fill = "red")+ #method = "lm"
  ylab("Mental health") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
  theme_bw()
Mental_Health = ggarrange(VT_descriptive,SF_descriptive,RE_descriptive,MH_descriptive)
ggsave(paste(savepath,"Mental_Health.pdf",sep=""),Mental_Health)
#
# check atrophy pattern by comparing the first (baseline) to final (10-years) timepoints
# tf0 = anat %>% filter(session == "0") %>% unique() %>% # select the right timepoints & unique data (no duplicates)
#   select(-eid, -age, -brainage, -corrected_brainage, -X, -session)
# tf120 = anat %>% filter(session == "120") %>% unique() %>% # select the right timepoints & unique data (no duplicates)
#   select(-eid, -age, -brainage, -corrected_brainage, -X, -session)
#
#
# CHECK TIME POINT DIFFERENCES (RATE OF CHANGE IN COHEN'S D)
# tf=anat %>% filter(session == "0" | session == "120") %>% unique() %>% # select the right timepoints & unique data (no duplicates)
#   select(-age, -brainage, -corrected_brainage, -X) %>% group_by(session) %>% filter( !duplicated(eid)) %>% group_split() #%>% select(eid, session) 
# tf=rbindlist(tf)
# a1 = data.frame(matrix(ncol=10,nrow=length(names(tf))))
# for (i in 1:length(names(tf))){
#   f1 = formula(paste(names(tf)[i],"~session",sep=""))
#   a1[i,] = tf %>% cohens_d(f1, paired = TRUE) %>%
#     inner_join(tf %>% t_test(f1, paired = TRUE, var.equal = FALSE)
#     )
# }
# names(a1) = names(tf %>% cohens_d(formula(paste(names(tf)[1],"~session",sep="")), paired = TRUE) %>%inner_join(tf %>% t_test(formula(paste(names(tf)[1],"~session",sep="")), paired = TRUE, var.equal = FALSE)))
# hist(a1$effsize,breaks = 20)
# hist(data$age,breaks=50)
#
#
#
#
#
#
#
# 1.2 Do any variables change significantly over time? ####
# Variables interesting for final analyses:
# Brain atrophy
# Brain age
# EDSS
# Cumulative relapses
# BL_Age
# BMI, BMI_change
# Omega3_suppl
# Treatment
# BL_BMI, rate of BMI change
# Vitamins
# lesion count & vol
# some data processing:
df = unique(data)
df = df %>% select(!starts_with("spm"))
# remove_c = c()
# for (i in 1:nrow(df)){
#   remove_c[i]=(ifelse(df$age[i+1] == df$age[i],"remove","keep"))
# }
#df$rem = remove_c
#df %>% filter(!remove_c == "remove") %>% filter(session == 120)
df1 = df %>% filter(session==0)
df2 = df %>% filter(session==120)
df1 = df1[!duplicated((df1)$eid),]
df2 = df2[!duplicated((df2)$eid),]
df1 = df1[df1$eid %in% df2$eid,]
df2 = df2[df2$eid %in% df1$eid,]
#df0 = rbind(df1,df2)
#
# test changes from baseline to 120 months
#t.test(df2$BAG_c,df1$BAG_c, paired = T) # sig increse in corrected and corrected BAG
#t.test(df2$BAG,df1$BAG, paired = T) # sig decrease in UNcorrected and uncorrected BAG

df2.1 = rbind(df1%>%select(eid,session,BAG_c),df2%>%select(eid,session,BAG_c))
df2.1$session = ifelse(df2.1$session == 0,0,120)
#cohens_d(df2.1,BAG_c~session,paired=T,ci=T)
df2.1 %>% group_by(session) %>% summarise(M = mean(na.omit(BAG_c)))

df1$BAG_c_diff = (df2$BAG_c-df1$BAG_c)
BAGdf = df1 %>% select(eid, BAG_c, BAG_c_diff)
#
demo10$edss_baseline = as.numeric(gsub(",",".",demo10$edss_baseline))
demo10$edss_month_12 = as.numeric(gsub(",",".",demo10$edss_month_12))
demo10$edss_month18 = as.numeric(gsub(",",".",demo10$edss_month18))
demo10$edss_month_24 = as.numeric(gsub(",",".",demo10$edss_month_24))
demo10$EDSSdiff_FUvsM24 = as.numeric(gsub(",",".",demo10$EDSSdiff_FUvsM24))
demo10$edss_month_120 = demo10$EDSSdiff_FUvsM24+demo10$edss_month_24
t.test(demo10$edss_month_120,demo10$edss_baseline,paired=T) # EDSS is higher at 120 months than at baseline
cohens_d(rbind(data.frame(session = "120", edss = demo10$edss_month_120),data.frame(session = "0", edss = demo10$edss_baseline)),edss~session,paired=T,ci=T,ref.group="120")
# showcasing the increase:
mean(na.omit(demo10$edss_month_120))
mean(na.omit(demo10$edss_baseline))
demo10$edss_ARC = (demo10$edss_month_120-demo10$edss_baseline)/(demo10$Age_OFAMS10-demo10$Age_BL_OFAMS)
# 
# Relapse
mean(na.omit(demo10$Cum_relapses))
t.test(na.omit(as.numeric(gsub(",",".",demo10$Relapserate_OFAMS))), mu = 0, alternative = "greater") # relapse rate > 0
t.test(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))), mu = 0, alternative = "greater") # total relapses > 0
min(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
max(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
median(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
mad(na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses))))
hist((na.omit(as.numeric(gsub(",",".",demo10$Cum_relapses)))),breaks = 20)
#
# Body mass index (BMI)
t.test(as.numeric(gsub(",",".",demo10$BL_BMI)), as.numeric(gsub(",",".",demo10$BMI_OFAMS10)), paired = T) # no sig differences in BMI
#
#
#
#
#
#
# EDSS plotting #######
demo10$eid=demo10$Patnr
edss_df0 = demo10%>% select(eid,Age_OFAMS10,edss_month_120,sex)
edss_df0$session = "144"
names(edss_df0) = c("eid", "age", "edss", "sex", "session")
edss_df = df3 %>% select(eid, age, edss, sex, session)
edss_df = rbind(edss_df,edss_df0)
# IMPORTANT: FILTER EMPTY DATA!!
edss_df = edss_df %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
edss_df = edss_df[order(edss_df$eid,edss_df$session),]
#write.csv(edss_df,paste(savepath,"tmp_table.csv",sep=""))
edss_df = read.csv(paste(savepath,"edss_age_table.csv",sep=""))
# edss_df=edss_df[order(edss_df$eid,as.numeri(edss_df$ession), edss_df$age),]
# df$age<-na.locf(df$age)
#edss_df %>% group_by(eid,session) %>% arrange(eid,session,age) %>% fill(age)
# 
# # we estimate missing ages
# a = (right_join(edss_df %>% filter(session == 0),demo10,by="eid")) %>% data.frame()
# a = ifelse((a)$age %>% is.na == F, a, (right_join(edss_df %>% filter(session == 6),demo10,by="eid")) %>% select(age) - 0.5)
# a = ifelse(a%>%is.na()== F, a, (edss_df %>% filter(session == 12))$age - 1)
# a = ifelse(a%>%is.na()== F, a, (edss_df %>% filter(session == 18))$age - 1.5)
# a = ifelse(a%>%is.na()== F, a, (edss_df %>% filter(session == 24))$age - 2)
# a = ifelse(a%>%is.na()== F, a, (edss_df %>% filter(session == 144))$age - 12)
# a1 = edss_df %>% filter(session == 0)
# a
# a1$age = a
# 
# b = ifelse((edss_df %>% filter(session == 6))$age %>% is.na == F, (edss_df %>% filter(session == 6))$age, (edss_df %>% filter(session == 0))$age + 0.5)
# b = ifelse(a%>%is.na()== F, b, (edss_df %>% filter(session == 12))$age - 0.5)
# b = ifelse(a%>%is.na()== F, b, (edss_df %>% filter(session == 18))$age - 1)
# b = ifelse(a%>%is.na()== F, b, (edss_df %>% filter(session == 24))$age - 1.5)
# b = ifelse(a%>%is.na()== F, b, (edss_df %>% filter(session == 144))$age - 11.5)
# b1 = edss_df %>% filter(session == 6)
# b1$age = b
# 
# c = ifelse((edss_df %>% filter(session == 12))$age %>% is.na == F, (edss_df %>% filter(session == 12))$age, (edss_df %>% filter(session == 0))$age + 1)
# c = ifelse(a%>%is.na()== F, c, (edss_df %>% filter(session == 6))$age + 0.5)
# c = ifelse(a%>%is.na()== F, c, (edss_df %>% filter(session == 18))$age - 0.5)
# c = ifelse(a%>%is.na()== F, c, (edss_df %>% filter(session == 24))$age - 1)
# c = ifelse(a%>%is.na()== F, c, (edss_df %>% filter(session == 144))$age - 11)
# c1 = edss_df %>% filter(session == 12)
# c1$age = c
# 
# d = ifelse((edss_df %>% filter(session == 18))$age %>% is.na == F, (edss_df %>% filter(session == 18))$age, (edss_df %>% filter(session == 0))$age + 1.5)
# d = ifelse(a%>%is.na()== F, d, (edss_df %>% filter(session == 6))$age + 1)
# d = ifelse(a%>%is.na()== F, d, (edss_df %>% filter(session == 12))$age + 0.5)
# d = ifelse(a%>%is.na()== F, d, (edss_df %>% filter(session == 24))$age - 0.5)
# d = ifelse(a%>%is.na()== F, d, (edss_df %>% filter(session == 144))$age - 10.5)
# d1 = edss_df %>% filter(session == 18)
# d1$age = d
# 
# e = ifelse((edss_df %>% filter(session == 24))$age %>% is.na == F, (edss_df %>% filter(session == 24))$age, (edss_df %>% filter(session == 0))$age + 2)
# e = ifelse(a%>%is.na()== F, e, (edss_df %>% filter(session == 6))$age + 1.5)
# e = ifelse(a%>%is.na()== F, e, (edss_df %>% filter(session == 12))$age + 1)
# e = ifelse(a%>%is.na()== F, e, (edss_df %>% filter(session == 18))$age + 0.5)
# e = ifelse(a%>%is.na()== F, e, (edss_df %>% filter(session == 144))$age - 10)
# e1 = edss_df %>% filter(session == 24)
# e1$age = e
# 
# f = ifelse((edss_df %>% filtfr(session == 144))$age %>% is.na == F, (edss_df %>% filter(session == 144))$age, (edss_df %>% filter(session == 0))$age + 12)
# f = ifelse(a%>%is.na()== F, f, (edss_df %>% filter(session == 6))$age + 11.5)
# f = ifelse(a%>%is.na()== F, f, (edss_df %>% filter(session == 12))$age + 11)
# f = ifelse(a%>%is.na()== F, f, (edss_df %>% filter(session == 18))$age + 10.5)
# f = ifelse(a%>%is.na()== F, f, (edss_df %>% filter(session == 24))$age + 10)
# f1 = edss_df %>% filter(session == 144)
# f1$age = f
# edss_df=rbind(a1,b1,c1,d1,e1,f1) #bind it
# edss_df = na.omit(edss_df)
# #
edss_df$session = factor(edss_df$session)
# to estimate the age and sex adjusted EDSS score, we use an LME
mml = lmer(edss ~ session + age + sex + (1|eid),edss_df[edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,])
# show the raw difference in EDSS score:
edss_df %>% filter(session == 144 & edss>=4) %>% summarize(mean(edss)) -edss_df[edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4))$eid,]%>% filter(session ==0) %>% summarize(mean(edss))
# get the time point 144 model coefficient
cofff=summary(mml)$coefficients[5]
# check also how many participants are left in stratified group
ids=edss_df[edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,]%>% select(eid) %>% unique  %>% unlist() %>% length
# plot it
e1 = ggplot(data= edss_df, aes(x = (session), y = edss)) + 
  geom_line(aes(group = eid), color = "grey50") + gghighlight(max(edss)>=4 ,use_direct_label = F) +
  geom_point(size = 1.5, alpha= 1, color = "grey50") + 
  stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("EDSS") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=3,y=8,label= paste("Adjusted 12-year difference =",round(cofff,2),"N =",ids),cex=4) + 
  ggtitle("Reached EDSS≥4")

mml2 = lmer(edss ~ session + age + sex + (1|eid),edss_df[!edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,])

#edss_df[edss_df$eid %in% (session ==144 & edss<4 | session ==24 & edss<4| session ==18 & edss<4| session ==12 & edss<4| session ==6 & edss<4)$eid,]
cofff=summary(mml2)$coefficients[5]
ids=(edss_df[!edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,])$eid%>%unique() %>% unlist() %>% length
level_order = c(0,6,12,18,24,144)
e2 = ggplot(data= edss_df, aes(x = (session), y = edss)) + 
  geom_line(aes(group = eid), color = "grey50") + gghighlight(max(edss)<4,use_direct_label = F) +
  geom_point(size = 1.5, alpha= 1, color = "grey50") + 
  stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("EDSS") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=3,y=8,label= paste("Adjusted 12-year difference =",round(cofff,2),"N =",ids),cex=4) + 
  ggtitle("Never reached EDSS≥4")
EDSS_descriptive = ggarrange(e1,e2)
# col = ifelse(edss_df$edss>=4,"red","gray85")
# EDSS_descriptive = ggplot(data= edss_df, aes(x = factor(session,level=level_order), y = edss)) + 
#   geom_point(size = 1.5, alpha= 1, color = "gray85") +
#   geom_path(aes(group = eid), color = "gray85") + #spaghetti plot
#   stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
#   ylab("EDSS") + xlab("Session (months)") +
#   stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.6)+ # geom = "pointrange"
#   theme_bw()
ggsave(paste(savepath,"EDSS_descriptive.pdf",sep=""),EDSS_descriptive, device = cairo_pdf,width = 9.5, height = 4)
#
#
# Specify characteristics of the functional decline group
ids = edss_df[edss_df$eid %in% (edss_df %>% filter(session ==144 & edss>=4 | session ==24 & edss>=4| session ==18 & edss>=4| session ==12 & edss>=4| session ==6 & edss>4))$eid,]%>% select(eid) %>% unique
edss_df$FLG = ifelse(edss_df$eid %in% unlist(ids), 1, 0)
#
#
#
# Paced auditory serial addition test (PASAT) assessment ####
t.test(demo10$PASAT_OFAMS10,demo10$BL_PASATcorrect,paired = T) # PASAT score decreases
# However, the simple t-test does not consider the fluctuations in PASAT scores over time - we need LMMs
pasat1 = demo10%>%select(Patnr,sex,BL_PASATcorrect,PASAT_24M, PASAT_OFAMS10)
pasat1 = melt(pasat1,id.vars = c("Patnr","sex"))
names(pasat1) = c("eid","sex","session","PASAT")
pasat1$session=as.numeric(pasat1$session)
pasat1$session=ifelse(as.numeric(pasat1$session) == 1,0,pasat1$session)
pasat1$session=ifelse(as.numeric(pasat1$session) == 2,24,pasat1$session)
pasat1$session=ifelse(as.numeric(pasat1$session) == 3,144,pasat1$session)
pasat1$session=factor(pasat1$session)
pasat1 = pasat1 %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
pasat1 = edss_df %>% right_join(pasat1, by = c("eid","session","sex"))

pasat1$FLG = ifelse(pasat1$eid %in% unlist(ids), 1, 0)
# 
# pmod = lmer(PASAT ~ session*FLG + age + sex + (1|eid),pasat1)
# summary(pmod)
# 


#FLG.PASAT = lm(PASAT ~ FLG + age + sex,pasat1%>%filter(session==0))
#summary(FLG.PASAT)
#
# Paced auditory serial addition test (PASAT) plotting ####
pasat1 %>% filter(session == 144 & PASAT<=40) %>% summarize(mean(PASAT)) -pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40))$eid,]%>% filter(session ==0) %>% summarize(mean(PASAT))
mml = lmer(PASAT ~ session + age + sex + (1|eid),pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40 | session ==24 & PASAT<=40| session ==0 & PASAT<=40))$eid,])
cofff=summary(mml)$coefficients[3] # get corrected coefficient / difference
# get ids 
pasat_ids=(pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=30 | session ==24 & PASAT<=30| session ==18 & PASAT<=30| session ==12 & PASAT<=30| session ==6 & edss>4))$eid,])$eid%>%unique() %>% unlist() %>% length
e1 = ggplot(data= pasat1, aes(x = factor(session,level=level_order), y = PASAT)) + 
  geom_line(aes(group = eid), color = "grey50") + gghighlight(min(PASAT)<=30,use_direct_label = F) +
  geom_point(size = 1.5, alpha= 1, color = "grey50") + 
  stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("PASAT") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=2,y=65,label= paste("Adjusted 10-year difference=",round(cofff,2),"N =",pasat_ids),cex=4) + 
  ggtitle("Reached PASAT≤40")
mml = lmer(PASAT ~ session + age + sex + (1|eid),pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT>40 | session ==24 & PASAT>40| session ==0 & PASAT>40))$eid,])
cofff=summary(mml)$coefficients[3]
pasat_ids=(pasat1[!pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=30 | session ==24 & PASAT<=30| session ==18 & PASAT<=30| session ==12 & PASAT<=30| session ==6 & edss>4))$eid,])$eid%>%unique() %>% unlist() %>% length
e2 = ggplot(data= pasat1, aes(x = factor(session,level=level_order), y = PASAT)) + 
  geom_line(aes(group = eid), color = "grey50") + gghighlight(min(PASAT)>30,use_direct_label = F) +
  geom_point(size = 1.5, alpha= 1, color = "grey50") + 
  stat_smooth(aes(group = 1),color="red",fill = "red", method = "lm")+ #
  ylab("PASAT") + xlab("Session (months)") +
  stat_summary(aes(group = 1), fun.data = "mean_cl_boot", shape = 17, size = 0.8,color="blue")+ # geom = "pointrange"
  theme_bw() + annotate("text",x=2,y=65,label= paste("Adjusted 10-year difference=",round(cofff,2),"N =",pasat_ids),cex=4) + 
  ggtitle("Reached PASAT≤40")

ggarrange(e1,e2)
# check N in stratified group
pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40 | session ==24 & PASAT<=40))$eid,] %>% group_by(session) %>% summarize(N = length(session))

#
pasat1[pasat1$eid %in% (pasat1 %>% filter(session ==144 & PASAT<=40 | session ==24 & PASAT<=40))$eid,] %>% select(eid) %>% unique
















# timed 25-foot walk test (T25FW) at 24 and 120 months
t.test(as.numeric(gsub(",",".",demo10$T_25FW_2)), as.numeric(gsub(",",".",demo10$T_25FW_1)), paired = T) # sig differences in T25FW
cohens_d(rbind(data.frame(session="24",T_25FW=as.numeric(gsub(",",".",demo10$T_25FW_1))),
               data.frame(session="120",T_25FW=as.numeric(gsub(",",".",demo10$T_25FW_2)))),
         T_25FW~session, ref.group="120", paired=T, ci=T)
demo10$T_25FW_ARC = (as.numeric(gsub(",",".",demo10$T_25FW_2))-as.numeric(gsub(",",".",demo10$T_25FW_1)))/14
#
# check proms (show no difference between BL and FU) [BL to 24 months]
df1 = df %>% filter(session==0)
df2 = df %>% filter(session==24)
df1 = df1[!duplicated((df1)$eid),]
df2 = df2[!duplicated((df2)$eid),]
df1 = df1[df1$eid %in% df2$eid,]
df2 = df2[df2$eid %in% df1$eid,]
t.test(df1$PF, df2$PF,paired=T)
t.test(df1$RF, df2$RF,paired=T)
t.test(df1$BP, df2$BP,paired=T)
t.test(df1$GH, df2$GH,paired=T)
t.test(df1$VT, df2$VT,paired=T)
t.test(df1$SF, df2$SF,paired=T)
t.test(df1$RE, df2$RE,paired=T)
t.test(df1$MH, df2$MH,paired=T)
#
# test total volume degeneration
icv %>% filter(session == 0 | session == 120) %>% t_test(TotalVol~session)
#icv %>% filter(session == 0 | session == 120) %>% cohens_d(TotalVol~session,ref.group="120", paired=T, ci=T)
t.test((icv %>% filter(session == 0))$TotalVol,(icv %>% filter(session == 120))$TotalVol)
#
vol.m = lmer(TotalVol~session+TIV+(1|eid),icv %>% filter(session == 0 | session == 120))
summary(vol.m) # estimate also the TIV corrected model (also reported in the article)
# df1 = anat %>% filter(session==0)
# df2 = anat %>% filter(session==120)
# df1 = df1[!duplicated((df1)$eid),]
# df2 = df2[!duplicated((df2)$eid),]
# df1 = df1[df1$eid %in% df2$eid,]
# df2 = df2[df2$eid %in% df1$eid,]
# df1$total.vol = rowSums(df1%>% select(contains("volume")))
# df2$total.vol = rowSums(df2 %>% select(contains("volume")))
# df1$vol_diff=(df2$total.vol-df1$total.vol)
# vols = df1 %>% select(eid,total.vol, vol_diff)
# t.test(rowSums(df2 %>% filter(session == 120) %>% select(contains("volume"))), rowSums(df1 %>% filter(session == 0) %>% select(contains("volume"))), paired=T) # clear decrease in brain vol
# df2=rbind(df1%>%select(session,total.vol),df2%>%select(session,total.vol))
# cohens_d(df2,total.vol~session,paired=T,ci=T)
#mean(rowSums(df1 %>% filter(session == 0) %>% select(contains("volume"))))
#mean((rowSums(df2 %>% filter(session == 120) %>% select(contains("volume")))))
#
# test fatigue changes
fati$fatigue = rowSums(fati[9:ncol(fati)])/9 # calculate fatigue sum score (higher = more)
fati$eid = fati$patno
fati1 = fati %>% filter(VISIT == "Month 24") %>% select(eid, fatigue)
fati2 = fati %>% filter(VISIT == "Baseline") %>% select(eid, fatigue)
fati1 = right_join(fati1,fati2,by="eid")
fati1$fatigue = fati1$fatigue.y-fati1$fatigue.x
names(fati1) = c("eid","fatigue_M24", "fatigue_BL", "fatigue_diff24")
#fati1 = fati1 %>% select(eid,fatigue)
t.test(fati1$fatigue_M24, fati1$fatigue_BL,mu = 0) # fatigue changes are not different from 0

# just double checking: relapses prior study (and meds) are not or even negatively related to relapse count in the followin 10 years
cor.test(demo10$relapses_12mnths_before_baseline,demo10$Cum_relapses,use="na.or.complete")
#
# correlations of edss and relapses over time [stable correlation between EDSS and relapses]
# relapses at the time reflect edss and vice versa at circa r = 0.20
# relapses 12 months before baseline 
# cor(demo10$relapses_12mnths_before_baseline, demo10$edss_baseline,use="na.or.complete")
# cor(demo10$relapses_12mnths_before_baseline, demo10$edss_month_12,use="na.or.complete")
# cor(demo10$relapses_12mnths_before_baseline, demo10$edss_month_18,use="na.or.complete")
# cor(demo10$relapses_12mnths_before_baseline, demo10$edss_month_24,use="na.or.complete")
# cor(demo10$relapses_12mnths_before_baseline, demo10$EDSS_10_short,use="na.or.complete")
#
# total/cumulative relapses
# cor(demo10$Cum_relapses, demo10$edss_baseline,use="na.or.complete")
# cor(demo10$Cum_relapses, demo10$edss_month_12,use="na.or.complete")
# cor(demo10$Cum_relapses, demo10$edss_month_18,use="na.or.complete")
# cor(demo10$Cum_relapses, demo10$edss_month_24,use="na.or.complete")
# cor(demo10$Cum_relapses, demo10$EDSS_10_short,use="na.or.complete")
#
# correlations of PROMS and relapses
prom$PF = prom%>%select(contains("spm3")) %>% rowSums() # Physical functioning (PF)
prom$RF = prom%>%select(contains("spm4")) %>% rowSums() # Role-physical (RF)
prom$BP = prom%>%select(spm7,spm8) %>% rowSums() # Bodily pain (BP)
prom$GH = prom%>%select(spm1,contains("spm11")) %>% rowSums() # General health (GH)
prom$VT = prom%>%select(spm9a,spm9b,spm9c,spm9d) %>% rowSums() # Vitality (VT)
prom$SF = prom%>%select(spm6,spm10) %>% rowSums() # Social functioning (SF)
prom$RE = prom%>%select(contains("spm5")) %>% rowSums() # Role-emotional (RE)
prom$MH = prom%>%select(spm9b, spm9c, spm9d, spm9f, spm3h) %>% rowSums() # Mental health (MH)
demo10$eid = demo10$Patnr
prom = prom %>% filter(session == 0) # select only baseline measure (there is anyways no change)
prom = right_join(demo10,prom, by="eid")
prom$BL_BMI = as.numeric(gsub(",",".",prom$BL_BMI))
prom$BMI_OFAMS10 = as.numeric(gsub(",",".",prom$BMI_OFAMS10))
prom$FatigueSS_score = as.numeric(gsub(",",".",prom$FatigueSS_score))
prom$Vit_A_0 = as.numeric(gsub(",",".",prom$Vit_A_0))
prom$Vit_D_0 = as.numeric(gsub(",",".",prom$Vit_D_0))
prom$Vit_E_0 = as.numeric(gsub(",",".",prom$Vit_E_0))
# hist(prom$Age_BL_OFAMS)
#
# lesion changes over time
# t.test(lesion_count$m120,lesion_count$baseline,paired=T)
# t.test(lesion_vol$m120,lesion_vol$baseline,paired=T)
# estimate cohen's d for changes
l = right_join(lesion_count %>% select(eid,baseline,m120),lesion_vol %>% select(eid,baseline,m120),by="eid")
names(l) = c("eid","baselineC","m120C","baselineV","m120V")
lesions = data.frame(LesionCount = c(l$m120C,l$baselineC), LesionVol=c(l$m120V,l$baselineV),
           Session = c(replicate(length(l$m120V),"120"),replicate(length(l$baselineV),"baseline")))
cohens_d(lesions,LesionCount~Session,ci = T,paired=T)
cohens_d(lesions,LesionVol~Session,ci = T,paired=T)
#
#
# blood changes over time
# NfL
t_test(formula=NfL..pg.ml.~Visit.nr,data = (blood%>%filter(Visit.nr==0 | Visit.nr == 24)),paired=T)
cohens_d(formula=NfL..pg.ml.~Visit.nr,data = (blood%>%filter(Visit.nr==0 | Visit.nr == 24)))
blood%>%filter(Visit.nr==0 | Visit.nr == 24)%>%group_by(Visit.nr)%>%summarize(M=mean(na.omit(NfL..pg.ml.)))
# CHI3L1
t_test(formula=CH3L.1..mg.ml..mean~Visit.nr,data = (blood%>%filter(Visit.nr==0 | Visit.nr == 24)),paired=T)
cohens_d(formula=CH3L.1..mg.ml..mean~Visit.nr,data = (blood%>%filter(Visit.nr==0 | Visit.nr == 24)))
blood%>%filter(Visit.nr==0 | Visit.nr == 24)%>%group_by(Visit.nr)%>%summarize(M=mean(na.omit(CH3L.1..mg.ml..mean)))
#
#
#
#
#remove duplicates
#full_join(icv %>% filter(session==0),icv %>% filter(session==120),by=c("eid","session"))
icv2 = icv %>% filter(session == "0" | session == "120") %>% unique()
remove_c = c()
for (i in 1:nrow(icv2)){ # loop provides Boolean-ish for exclusion based on duplicates (seciond is being kept)
  remove_c[i]=(ifelse(icv2$eid[i] == icv2$eid[i-1],F,T)) # T/F = keep/remove
}
remove_c = ifelse(is.na(remove_c) == T, T,remove_c) # one row in beginning must be NA due to loop implementation
icv = icv2[remove_c,]
icv2 = icv%>%filter(session==0)
icv3 = icv%>%filter(session==120)
icv_change = right_join(icv2,icv3,by="eid")
icv_change$TotalVol_change=(icv_change$TotalVol.y-icv_change$TotalVol.x)
icv_change$TIV_change=(icv_change$TIV.y-icv_change$TIV.x)
#
# prom = right_join(vols,prom,by="eid") # proms
# prom$vol_ARC = prom$vol_diff/(prom$Age_OFAMS10-prom$Age_BL_OFAMS)
prom = right_join(BAGdf,prom,by="eid")
#prom$BAG_ARC = prom$BAG_c_diff/(prom$Age_OFAMS10-prom$Age_BL_OFAMS)
# mat1 = prom %>% select(Age_BL_OFAMS,BAG_c, BAG_ARC, edss_ARC, pasat_ARC, PASAT_24M, 
#                        edss_baseline, edss_month_12, edss_month18, edss_month_24, edss_month_120,
#                        relapses_12mnths_before_baseline,
#                        Cum_relapses, PF,RF,BP,GH,VT,SF,RE,MH,
#                        BL_BMI,FatigueSS_score,mean_vitA, mean_vitD,mean_vitE
# ) %>% cor(use="na.or.complete")
# #Omega3_suppl, Treatment_OFAMS, Smoke_last_10y
# pdf(paste(savepath,"cormat.pdf",sep=""), width=18, height=18)
# heatmap.2(Rowv = F, x = mat1, col = col, symm = T, cellnote=round(mat1,2),notecol="black", dendrogram = "none",key=T,density.info="density",trace = "none",notecex = 1.25)
#heatmap.2(Rowv = F, x = mat1, col = col, symm = T, cellnote=round(mat1,2),notecol="black",notecex = 1.25)
# dev.off()
#
#
# Nb of relapses is weakly reflected by proms
# Nb of relapses is acceptably (but not strongly) reflected by EDSS when later in time
# Previous relapses indicate further relapses in the future
# Age associates with everything
# BAGc

# plot(prom$BAG_c,prom$Age_BL_OFAMS)
# plot(prom$edss_ARC,prom$Age_BL_OFAMS)
# plot(prom$edss_ARC,prom$Age_BL_OFAMS)
# plot(prom$edss_ARC,prom$BAG_ARC) # 
# cor(prom$BAG_c,(prom$pasat_ARC),use="na.or.complete")
# cor(prom$BAG_c,(prom$edss_ARC),use="na.or.complete")

# mean(na.omit(demo10$edss_baseline))
# mean(na.omit(demo10$edss_month_120))
#
#
# Analyses across variables (simple models) ####
# 2. Explain annual rate of change (ARC) ####
prom=right_join(prom,geno,by="eid") # merge relevant dfs
# prom$Relapserate_OFAMS = as.numeric(gsub(",",".",prom$Relapserate_OFAMS)) # add relapse rate
prom$Mean_BMI_OFAMS = as.numeric(gsub(",",".",prom$Mean_BMI_OFAMS)) # add bmi
prom$age_gap = prom$Age_OFAMS10-prom$Age_BL_OFAMS # add age gap
prom = right_join(icv2,prom,by="eid")
prom = right_join(l,prom,by="eid")
BL_blood = blood %>% filter(Visit.nr==0)
BL_blood$eid=BL_blood$Sub
prom=right_join(prom,BL_blood,by="eid")

# 2.1 All Ordinary least squares ####
# Reasoning: identify most influential variables
# For that, use leave one out CV
#
#
############ PASAT
# select variables for first model
# pasat_df = prom%>%dplyr::select(pasat_ARC,sex, geno, Age_BL_OFAMS, 
#                          relapses_12mnths_before_baseline,edss_baseline,
#                          smoking_OFAMS,BL_BMI,BAG_c,TotalVol,TIV,
#                          baselineC,baselineV,
#                          ,PF,RF,BP,GH,VT,SF,RE,MH)
# # impute missing values
# pasat_df = impute_shd(pasat_df, .~1, backend="VIM")
# #apply(Y, c(1, 2), mean, na.rm = TRUE)
# a = data.frame(Counter = 1:length(cvfit$cvm),Lambda = cvfit$lambda, Error = cvfit$cvm) # make a df of the obtained models' errors and lambdas
# a = a[order(a$Error),][1:10,] # restrict to the best 20 cross-validated models
# a = rowSums(data.frame(cvfit$coeff)[,a$Counter])/10 # average across the top 20
# a = data.frame(Name=names(a),Weights=as.numeric(a))%>%filter(abs(Weights)>0.009)
# f1 = paste("pasat_ARC~",paste(a$Name[2:length(a$Name)],collapse ="+")) # make into formula
# mod1 = lm(f1,data=pasat_df)
# summary(mod1)
#
#
#
############ PASAT
# select variables for first model
# pasat_df = prom%>%dplyr::select(edss_ARC,sex, geno, Age_BL_OFAMS, 
#                                 relapses_12mnths_before_baseline,PASAT_24M,
#                                 smoking_OFAMS,BL_BMI,BAG_c,TotalVol,TIV,
#                                 baselineC,baselineV,
#                                 ,PF,RF,BP,GH,VT,SF,RE,MH)
# impute missing values
# pasat_df = impute_shd(pasat_df, .~1, backend="VIM")
#
#prom$FLG = ifelse(prom$eid %in% unlist(ids), 1, 0)

prom=prom%>%select(-Sex_OFAMS10,-sex)
jj=left_join(edss_df %>% filter(session==0),pasat1 %>% filter(session==0) %>% select(eid,PASAT), by="eid")
prom=left_join(jj,prom,by="eid")
prom = prom %>% select(eid,FLG,age,sex,edss,PASAT,geno,relapses_12mnths_before_baseline,
                       CH3L.1..mg.ml..mean,NfL..pg.ml.,smoking_OFAMS,Mean_BMI_OFAMS,
                       Treatment_OFAMS,Omega3_suppl,baselineC,baselineV,
                       PF,RF,BP,GH,VT,SF,RE,MH,Vit_A_0,Vit_D_0,Vit_E_0)


prom=left_join(df2 %>% select(eid,BAG_c),prom,by="eid")
write.csv(prom,paste(savepath,"interrim_data.csv",sep=""))
# run all models
# m1=lm(pasat_ARC~sex+geno+Age_BL_OFAMS +relapses_12mnths_before_baseline+
#         edss_baseline +
#         smoking_OFAMS+Mean_BMI_OFAMS+
#         Treatment_OFAMS+BAG_c+TotalVol+TIV+
#         baselineC + baselineV +
#         +PF+RF+BP+GH+VT+SF+RE+MH, data=prom)
# #numCores = detectCores() 
# pasat_all = ols_step_all_possible(m1,max_order = 10) #
# pasat_all
# plot(pasat_all) # plot the different models
#
#
############### EDSS
# m1=lm(edss_ARC~sex+geno+Age_BL_OFAMS +pasat_ARC+PASAT_24M+
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+BAG_c+TotalVol+TIV+T_25FW_1+
#         baselineC + baselineV +
#         +PF+RF+BP+GH+VT+SF+RE+MH+relapses_12mnths_before_baseline, data=prom)
# edss_all = ols_step_all_possible(m1)
# pasat_all
# plot(pasat_all)
# 
# # 2.2 Linear Models ####
# #
# m1=lm(pasat_ARC~sex+geno+Age_BL_OFAMS +relapses_12mnths_before_baseline+
#         edss_baseline +
#         smoking_OFAMS+Mean_BMI_OFAMS+
#         Treatment_OFAMS+BAG_c+TotalVol+TIV+
#         baselineC + baselineV +
#         +PF+RF+BP+GH+VT+SF+RE+MH, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# # m1=lm(pasat_ARC~sex+geno+Age_BL_OFAMS + age_gap +relapses_12mnths_before_baseline+
# #         smoking_OFAMS+Mean_BMI_OFAMS+Relapserate_OFAMS+
# #         Treatment_OFAMS+Omega3_suppl+BAG_c+vol_ARC+
# #         +PF+BP+GH+VT+SF, data=prom)
# #
# # winning model explaining PASAT
# m2=lm(pasat_ARC~geno+Age_BL_OFAMS+BP+RE+baselineC,data=prom)
# vif(m2,type = "predictor")
# vif(m2,type = "terms")
# summary(m2)
# model_parameters(m2, bootstrap = TRUE, iterations = 10000)
# #
# #
# prom$T_25FW_1=as.numeric(gsub(",",".",prom$T_25FW_1))
# # first iteration explaining future edss changes
# m1=lm(edss_ARC~sex+geno+Age_BL_OFAMS +pasat_ARC+PASAT_24M+
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+BAG_c+TotalVol+TIV+T_25FW_1+
#         baselineC + baselineV +
#         +PF+RF+BP+GH+VT+SF+RE+MH+relapses_12mnths_before_baseline, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# # second
# m1=lm(edss_ARC~Treatment_OFAMS+
#         Omega3_suppl+TotalVol+TIV+T_25FW_1+
#         +PF+RF+VT+SF+RE+MH, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# #
# #
# 
# m1=lm(edss_ARC~sex+Age_BL_OFAMS +
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+
#         +PF+RF+VT+SF+RE+MH+relapses_12mnths_before_baseline, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# 
# m1=lm(edss_ARC~sex+Age_BL_OFAMS +
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+
#         +PF+RF+BP+VT+SF+RE+relapses_12mnths_before_baseline, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# m1=lm(edss_ARC~sex+Age_BL_OFAMS +
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+
#         +PF+RF+BP+VT+SF+RE, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# m1=lm(edss_ARC~sex+Age_BL_OFAMS +
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+
#         +PF+RF+VT+SF+RE, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# m1=lm(edss_ARC~Age_BL_OFAMS +
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+
#         +PF+RF+VT+SF+RE, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# # the winning model for edss
# m1=lm(edss_ARC~
#         smoking_OFAMS+ BL_BMI+Treatment_OFAMS+
#         Omega3_suppl+
#         +PF+RF+VT+SF+RE, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# model_parameters(m1, bootstrap = TRUE, iterations = 10000)
# #
# # difficult to explain relapses
# m1=lm(Cum_relapses~sex+geno+Age_BL_OFAMS + PASAT_24M+
#         age_gap +smoking_OFAMS+ Mean_BMI_OFAMS +
#         Treatment_OFAMS+Omega3_suppl+BAG_c+TotalVol+TIV+
#         edss_baseline+relapses_12mnths_before_baseline,
#       data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# m1=lm(Cum_relapses~geno+Age_BL_OFAMS + 
#         age_gap +smoking_OFAMS+ Mean_BMI_OFAMS +
#         Omega3_suppl+TotalVol+TIV,
#       data=prom)
# summary(m1)
# m1=lm(Cum_relapses~geno+Age_BL_OFAMS + 
#         smoking_OFAMS+ Mean_BMI_OFAMS +
#         Omega3_suppl+TotalVol+TIV,
#       data=prom)
# summary(m1)
# 
# m1=lm(Cum_relapses~geno + 
#         smoking_OFAMS+ Mean_BMI_OFAMS +
#         Omega3_suppl,
#       data=prom)
# summary(m1)
# 
# # as well as the rate of relapses
# m1=lm(Relapserate_OFAMS~sex+geno+Age_BL_OFAMS + age_gap +
#         smoking_OFAMS+ Mean_BMI_OFAMS +Treatment_OFAMS+
#         Omega3_suppl+BAG_c+TotalVol+TIV+
#         edss_baseline+relapses_12mnths_before_baseline, data=prom)
# summary(m1)
# m1=lm(Relapserate_OFAMS~
#         smoking_OFAMS+ Mean_BMI_OFAMS +Treatment_OFAMS+
#         Omega3_suppl+BAG_c+TotalVol+TIV+
#         edss_baseline+relapses_12mnths_before_baseline, data=prom)
# summary(m1)
# #
# #
# # finally walking test
# m1=lm(T_25FW_ARC~sex+geno+Age_BL_OFAMS +relapses_12mnths_before_baseline+
#         edss_baseline +PASAT_24M+
#         smoking_OFAMS+Mean_BMI_OFAMS+
#         Treatment_OFAMS+Omega3_suppl+BAG_c+TotalVol+TIV+
#         +PF+RF+BP+GH+VT+SF+RE+MH, data=prom)
# vif(m1,type = "predictor")
# vif(m1,type = "terms")
# summary(m1)
# m1=lm(T_25FW_ARC~sex+relapses_12mnths_before_baseline+
#         edss_baseline +
#         smoking_OFAMS+
#         BAG_c+
#         +PF+GH, data=prom)
# summary(m1)
# m1=lm(T_25FW_ARC~edss_baseline +
#         smoking_OFAMS+
#         BAG_c+
#         +PF+GH, data=prom)
# summary(m1)
#
# #Some test on outliers
# prom12 = prom%>%filter(Cum_relapses>1)
# plot(prom12$Age_OFAMS10,prom12$vol_ARC)
# 
# plot(prom$edss_ARC, prom$vol_ARC)
# plot(prom12$edss_ARC, prom12$vol_ARC)
# 
# cor(prom$edss_ARC, prom$vol_ARC,use="na.or.complete")
# cor(prom12$edss_ARC, prom12$vol_ARC,use="na.or.complete")
# 
# exp_m1=lm(edss_ARC~vol_ARC*Cum_relapses,data=prom)
# summary(exp_m1)
#
#
#
# 3. Longitudinal analyses (mixed models) ####
########### CHECK LONGITUDINAL CHANGES FOR BRAIN AGE AND EDSS
#
#
# modelling plans:
# a) check what can explain BAG (changes)
# b) check the influence of relapses and other statics
# All of this expressed in formulas:
# a) lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+GH+VT+SF+RE+MH+(1|eid),data=data)
# b) lm(BAG_c_change~edss_change+sex+prom_change+nb_of_relapses, data = lm_dat)
# for the final model: add "DISEASE_DURATION"   "SYMPTOM_DURATION"   "AGE_SYMPTOM_ONSET"  "AGE_DIAGNOSIS"?
# also: onset or final age? all of those are in the demo data
#
#
#
# 3.1 Brain age changes ####
# mod1 = lmer(BAG_c~session+sex+geno+age+(age||eid),data=data) # Convergence issues
# mod = lmer(BAG_c~session+sex+geno+age+(1|eid)+(0+age|eid),data=data) # Convergence issues
mod1 = lmer(BAG_c~session+sex+geno+age+(0+age|eid),data=data)
mod2 = lmer(BAG_c~session+sex+geno+(0+age|eid),data=data)
mod3 = lmer(BAG_c~session+sex+geno+age+(1|eid),data=data)
anova(mod1,mod2,mod3) # mod3 is better than the other models
#
summary(mod3) # geno and age seem interesting
ggpredict(mod3, terms=c("sex"), type = "fe")
ggpredict(mod3, terms=c("geno"), type = "fe")
ggpredict(mod3, terms=c("age"), type = "fe")
ggpredict(mod3, terms=c("sex","age"), type = "fe")
ggpredict(mod3, terms=c("sex","geno"), type = "fe") 
ggpredict(mod3, terms=c("age","geno"), type = "fe") # geno seems interesting, but likely no interaction with age
#
ri_m_plot = sjPlot::plot_model(mod3, pred.type="re", ci.lvl=.95)
ri_m_plot = ri_m_plot +ggtitle("Corrected BAG Predictors")+ theme_bw()
ggsave(paste(savepath,"BAG_Predictors.pdf",sep=""),ri_m_plot)
#
sjPlot::plot_model(mod3, type="pred", terms=c("session","geno"), pred.type="fe", ci.lvl=.95)
sjPlot::plot_model(mod3, type="pred", terms=c("age","geno"), pred.type="fe", ci.lvl=.95)

# add clinical vars for explanations
mod4 = lmer(BAG_c~session+sex+geno+age+edss+(1|eid),data=data)
summary(mod4) # edss is little indicative of BAG
# sjPlot::plot_model(mod4, pred.type="re", ci.lvl=.95)
# sjPlot::plot_model(mod4, type="pred", terms=c("age"), pred.type="fe", ci.lvl=.99)
# sjPlot::plot_model(mod4, type="pred", terms=c("session"), pred.type="fe", ci.lvl=.99)
# sjPlot::plot_model(mod4, type="pred", terms=c("session","geno"), pred.type="fe", ci.lvl=.99)
#
# check for PROMs
mod5.1 = lmer(BAG_c~session+sex+geno+age+PF+(1|eid),data=data)
mod5.2 = lmer(BAG_c~session+sex+geno+age+edss+PF+(1|eid),data=data)
mod6 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+(1|eid),data=data)
mod7 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+(1|eid),data=data) # only proper improvement
mod8 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+GH+(1|eid),data=data)
mod9 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+GH+VT+(1|eid),data=data)
mod10 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+GH+VT+SF+(1|eid),data=data)
mod11 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+GH+VT+SF+RE+(1|eid),data=data)
mod12 = lmer(BAG_c~session+sex+geno+age+edss+PF+RF+BP+GH+VT+SF+RE+MH+(1|eid),data=data)
anova(mod5.1,mod5.2,mod6,mod7,mod8,mod9,mod10,mod11,mod12)
summary(mod7)
sjPlot::plot_model(mod7, pred.type="re", ci.lvl=.95) # BP effect the wrong way, false positive
# After all, only genotype matters! (no proms, no edss)
#
# 3.2 Brain change ####
# 3.2.1 Prep ####
# percentage in brain change
tf=anat %>% filter(session == "0" | session == "120") %>% unique() %>% # select the right timepoints & unique data (no duplicates)
  select(-age, -brainage, -corrected_brainage, -X) %>% group_by(session) %>% filter( !duplicated(eid)) %>% group_split()
tf1 = tf[[1]]
tf2 = tf[[2]]
tf1 = tf1[tf1$eid %in% tf2$eid,]
tf2 = tf2[tf2$eid %in% tf1$eid,]
tf1 = tf1[order(tf1$eid),] # order data
tf2 = tf2[order(tf2$eid),]
eidlist = tf1$eid
tfdf=(tf2-tf1)/tf1
tfdf$eid = eidlist
edss$change = (edss$M24_EDSSscore-edss$BL_EDSSscore) # estimate absolute edss change in the first 24 months
demo$eid = demo$Patnr
# merge all data
tfdf = right_join(tfdf,geno,by="eid")
tfdf = right_join(tfdf,fati1,by="eid")
tfdf = right_join(tfdf,demo,by="eid")
#
# get more data for modelling
# moredat = na.omit(data) %>% filter(session == 0 | session == 24) %>% unique()
# more1 = (moredat[moredat$eid %in% eidlist,] %>% filter(session == 24) %>% select(eid,PF,RF,BP,GH,VT,SF,RE,MH))
# more2 = (moredat[moredat$eid %in% more1$eid,] %>% filter(session == 0) %>% select(eid,PF,RF,BP,GH,VT,SF,RE,MH))
# more1 = more1[more1$eid %in% more2$eid,]
# more = more2-more1
# more$eid = more1$eid
# tfdf = merge(tfdf, more, by = "eid")
#
#
#
# 3.2.2 Plot raw percentage change 120 months after baseline ####
plotprep = function(Volume_Surf_Thick){
  vols = tfdf %>% select(contains(Volume_Surf_Thick))
  vals = c()
  for (i in 1:ncol(vols)){
    vals[i] = mean(na.omit(vols[,i]))
  }
  ggdf = data.frame(RegNames = names(vols), Change = vals)
  ggdf = ggdf[order(ggdf$RegNames),]
  ggdf$hemi = ifelse(grepl("lh_",ggdf$RegNames)==T,"left","right")
  #
  brain_labels(dk)[!brain_labels(dk) %in% gsub("_volume","",ggdf$RegNames)] # CC not included in our pipeline. Removed for plotting
  # Hence: this has to be fixed afterwards
  return(ggdf)
}
ggvol = plotprep("volume")
ggvol$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",ggvol$RegNames)][1:34]
ggsurf = plotprep("area")
ggsurf$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_area","",ggsurf$RegNames)][1:35]
ggthick = plotprep("thick")
ggthick$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_thickness","",ggthick$RegNames)][1:35]
p1 = ggplot(ggvol) + geom_brain(atlas = dk,aes(fill = Change),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  labs(title="Raw percentual volume changes over 120 months") + 
  theme_void()
p2 = ggplot(ggsurf) + geom_brain(atlas = dk,aes(fill = Change),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") + 
  labs(title="Raw percentual surface area changes over 120 months") + 
  theme_void()
p3 = ggplot(ggthick) + geom_brain(atlas = dk,aes(fill = Change),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") + #"firebrick","goldenrod"
  labs(title="Raw percentual thickness changes over 120 months") + 
  theme_void()
set_scale_union(p1,p2,p3,scale=scale_fill_gradient2(low = "blue",mid = "white",high="red"))
raw_percentage_change = ggarrange(p1,p2,p3, ncol=1, common.legend = T, legend = "bottom")
#
# 
# 3.2.3 Plot effect sizes for adjusted change 120 months after baseline ####
# brain change adjusted by sex, age, geno, Treatment, BMI
# for that, merge demo and anat data frames
dmg = demo10 %>% select(eid, BL_BMI, BMI_OFAMS10) # get BMI
dmg$BL_BMI = as.numeric(gsub(",",".",dmg$BL_BMI))
dmg$BMI_OFAMS10 = as.numeric(gsub(",",".",dmg$BMI_OFAMS10))
dmg = melt(dmg,id.vars = "eid")
dmg$variable = ifelse(dmg$variable == "BL_BMI",0,120)
names(dmg) = c("eid","session","BMI")
anat = full_join(dmg,anat,by=c("eid","session"))
# add smoking (as constant)
demo10$smoked = ifelse(ifelse(is.na(demo10$Smoke_last_10y)==T,0,demo10$Smoke_last_10y) + ifelse(is.na(demo10$smoking_OFAMS)==T,0,demo10$smoking_OFAMS)>0,1,0)
anat = right_join(anat,demo10%>%select(eid,smoked),by="eid")
# add constants
dmg = demo %>% select(eid, HLA_1501_1, Gender,Treatment)
dmg$session = 0
dmg2 = dmg
anat = right_join(dmg %>% select(eid, HLA_1501_1),anat,by=c("eid"))
anat = right_join(dmg %>% select(eid, Gender),anat,by=c("eid"))
anat = right_join(dmg %>% select(eid, Treatment),anat,by=c("eid"))
anat$session=as.factor(anat$session)
anat = full_join(anat,data%>%select(eid,session,edss),by=c("eid","session"))
anat = (unique(anat))
anat = anat %>%group_by(eid) %>%fill(Gender,HLA_1501_1,Treatment)
imp = impute_shd(anat, edss~1, backend="VIM") # hence imp is left out
imp$BAG = imp$corrected_brainage-imp$age
imp$age = scale(imp$age)
imp$BMI = scale(imp$BMI)
imp$eid = factor(imp$eid)
metrics = c("volume","area","thickness")
imp = imp%>%select(!contains("Mean")) %>% select(!contains("WhiteSurfArea"))
regs=mods=list()
for (metric in 1:length(metrics)){
  brain_vars=imp %>% dplyr::select(starts_with("lh_")|starts_with("rh_"))%>% select(contains(metrics[metric])) %>%names()
  #brain_vars = brain_vars[1:length(brain_vars)]
  age=sex=gene=bmi=treat=smoke=BAG=edss=r2m=r2c=c()
  p=data.frame(matrix(nrow=length(brain_vars),ncol=10)) # output data frame for p-vals: ncol=8 due to nb of vars of interest
  for (i in 1:length(brain_vars)){
    imp[brain_vars][i] = scale(imp[brain_vars][i])
    f1 = formula(paste(brain_vars[i],"~age+Gender+HLA_1501_1+BMI+
                       Treatment+smoked+BAG+edss+(1|eid)"))
    model = lmer(f1,imp)
    age[i]=summary(model)$coefficients[2]
    sex[i]=summary(model)$coefficients[3]
    gene[i]=summary(model)$coefficients[4]
    bmi[i]=summary(model)$coefficients[5]
    treat[i]=summary(model)$coefficients[6]
    smoke[i]=summary(model)$coefficients[7]
    BAG[i]=summary(model)$coefficients[8]
    edss[i]=summary(model)$coefficients[9]
    r2m[i] = r.squaredGLMM(model)[1]
    r2c[i] = r.squaredGLMM(model)[2]
    p[i,] = c(brain_vars[i],summary(model)$coefficients[,5])
    names(p) = c("Region","Intercept","Age","Sex","HLA_1501","BMI","Treatment","Smoked","BAG","EDSS")
    # reg[i]=(predictions(model))$estimate # note, this is ignoring the random effect!
    #for assessing the average prediction, one can also use avg_predictions(model)
  }
  regs[[metric]] = data.frame(RegNames = brain_vars, Age = age, 
                              Sex = sex, HLA15011 = gene, 
                              BMI = bmi, Treatment = treat,
                              Smoking = smoke, BAG, EDSS = edss,
                              R2m=r2m,R2c=r2c)
  mods[[metric]] = p
}
multiplot=function(effect){
  for (i in 1:length(regs)){
    regs[[i]]$hemi = ifelse(grepl("lh_",regs[[i]]$RegNames)==T,"left","right")
  }
  #brain_labels(dk)[!brain_labels(dk) %in% gsub("_volume","",regs[[1]]$RegNames)] # CC not included in our pipeline. Removed for plotting
  regs[[1]]$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",regs[[1]]$RegNames)][1:34]
  regs[[2]]$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_area","",regs[[2]]$RegNames)][1:34]
  regs[[3]]$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_thickness","",regs[[3]]$RegNames)][1:34]
  p1=ggplot(regs[[1]] %>% select(region,hemi,RegNames, !!!syms(effect)))+
    geom_brain(atlas = dk,aes_string(fill = paste(effect)),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Effect of ",paste(effect)," on Volume",sep="")) + 
    theme_void() + labs(fill="Effect size")
  p2=ggplot(regs[[2]] %>% select(region,hemi,RegNames, !!!syms(effect)))+
    geom_brain(atlas = dk,aes_string(fill = paste(effect)),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Effect of ",paste(effect)," on Cortical Surface Area",sep="")) + 
    theme_void() + labs(fill="Effect size")
  p3=ggplot(regs[[3]] %>% select(region,hemi,RegNames, !!!syms(effect)))+
    geom_brain(atlas = dk,aes_string(fill = paste(effect)),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Effect of ",paste(effect)," on Cortical Thickness",sep="")) + 
    theme_void() + labs(fill="Effect size")
  set_scale_union(p1,p2,p3,scale=scale_fill_gradient2(low = "blue",mid = "white",high="red"))
  return(ggarrange(p1,p2,p3, ncol=1, common.legend = T, legend = "bottom"))
}
plotlist=list()
for (i in 1:length(names(regs[[metric]])[2:9])){
  plotlist[[i]]=multiplot(names(regs[[metric]])[2:9][i])
}
# just as a reminder: sex=female, HLA150=present,smoking=yes
# plot all
ggarrange(plotlist = plotlist)
# plot adjusted p-vals
nbtest=length(metrics)*length(brain_vars) # total number of tests = alpha adjustment
0.05/nbtest # adjusted alpha
vol.padj=area.padj=thick.padj=mods[[1]]
for (i in 1:length(mods[[1]])-1){
  vol.padj[i+1] = as.numeric(unlist(mods[[1]][i+1]))*nbtest
  area.padj[i+1] = as.numeric(unlist(mods[[2]][i+1]))*nbtest
  thick.padj[i+1] = as.numeric(unlist(mods[[3]][i+1]))*nbtest
}
vol.padj$Region=area.padj$Region=thick.padj$Region=mods[[1]]$Region
#beta_names=names(vol.padj)[3:ncol(vol.padj)]

# to do: dense down data frame to contain only variables of interest
# then, write all p-vals > 1 as 1
# then, plot the p-vals per beta, per brain feature (surf,area,thick)
ifelse(vol.padj[3:ncol(vol.padj)]>1,1,vol.padj[3:ncol(vol.padj)])

vol.padj[names(vol.padj)[3:ncol(vol.padj)]] %>% filter_all(any_vars(. < .05))
#
# Discontinued and other snippets ####
# Perform repeated cross-validation to find optimal lambda value
# set.seed(1234)
# ctr = trainControl(method = "repeatedcv", # define CV procedure (10 folds, 1000 reps)
#                    repeats = 1000,
#                    verboseIter = TRUE
#                    #seeds = seeds
#                    )
# x = pasat_df %>% dplyr::select(-pasat_ARC) %>% as.matrix()
# y = pasat_df %>% dplyr::select(pasat_ARC) %>% unlist() %>% as.numeric
# # run generalized linear model using elastic net penalty
# cvfit = cvr.glmnet(X=x,Y=y,family="gaussian",standardize=T,nfolds=3,ncv=1000,type.measure="mse")
# impute missing values
#
#
#
#
#
# pasat_df = impute_shd(pasat_df, .~1, backend="VIM")
# x = pasat_df %>% dplyr::select(-edss_ARC) %>% as.matrix()
# y = pasat_df %>% dplyr::select(edss_ARC) %>% unlist() %>% as.numeric
# # run generalized linear model using elastic net penalty
# cvfit = cvr.glmnet(X=x,Y=y,family="gaussian",standardize=T,nfolds=3,ncv=1000,type.measure="mse")
# #apply(Y, c(1, 2), mean, na.rm = TRUE)
# a = data.frame(Counter = 1:length(cvfit$cvm),Lambda = cvfit$lambda, Error = cvfit$cvm) # make a df of the obtained models' errors and lambdas
# a = a[order(a$Error),][1:10,] # restrict to the best 20 cross-validated models
# a = rowSums(data.frame(cvfit$coeff)[,a$Counter])/10 # average across the top 20
# a = data.frame(Name=names(a),Weights=as.numeric(a))%>%filter(abs(Weights)>0.009)
# f1 = paste("edss_ARC~",paste(a$Name[2:length(a$Name)],collapse ="+")) # make into formula
# mod1 = lm(f1,data=pasat_df)
# summary(mod1)