# Explainability / Feature importance for model used on OFAMS data
#
# Max Korbmacher, 30, Nov 2024
#
#
# Packages and data ####
#
# pkgs
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,ggseg,dplyr,ggpubr)
#
# data
# (note that all the data used here can be found in the associated OSF repo: https://osf.io/3f4md/)
setwd("/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/Hemi2/Explain_model")
filenames = list.files(path="explain/",pattern="*0_BOTH.csv", full.names=TRUE) # stratified knockout importance for both-hemisphere model (knockout_res*.csv)
ldf = lapply(filenames, read.csv)
filenames = list.files(path="explain/",pattern="*0_RIGHT.csv", full.names=TRUE) # stratified knockout importance for both-hemisphere model (knockout_res*.csv)
ldfR = lapply(filenames, read.csv)
filenames = list.files(path="explain/",pattern="*0_LEFT.csv", full.names=TRUE) # stratified knockout importance for both-hemisphere model (knockout_res*.csv)
ldfL = lapply(filenames, read.csv)
filenames = list.files(path="FeatureImp",pattern="*.csv", full.names=TRUE) # feature weights for each model (weights*.csv)
fdf = lapply(filenames, read.csv)
# filenames = list.files(path="explain",pattern="*.csv", full.names=TRUE) #  knockout importance
# explain = lapply(filenames, read.csv)
kB=read.csv("explain/knockout_BOTH.csv")
explain=list(kB)
#
# Functions ####
# Make plot scale unifier function
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
# Make a plotting functions
## for feature weight based importance
impp=function(data,Hemisphere){
  data$hemi=ifelse(grepl("lh_",data$Parameter)==T,"left","right")
  data = data[order(data$Parameter),]
  vols=data[grepl("_volume",data$Parameter),]
  thick=data[grepl("_thickness",data$Parameter),]
  surf=data[grepl("_area",data$Parameter),]
  vols$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(gsub("_volume","",vols$Parameter),3)][1:34]
  thick$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(gsub("_thickness","",thick$Parameter),3)][1:34]
  surf$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(gsub("_area","",surf$Parameter),3)][1:34]
  #
  plot1=ggplot(vols) + geom_brain(atlas = dk,aes(fill = F),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Feature Importance of Volumes", sep="")) + labs(fill='F') +
    theme_void()
  plot2=ggplot(thick) + geom_brain(atlas = dk,aes(fill = F),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Feature Importance of Cortical Thickness", sep="")) + labs(fill='F') +
    theme_void()
  plot3=ggplot(surf) + geom_brain(atlas = dk,aes(fill = F),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Feature Importance of Surface Area", sep="")) + labs(fill='F') +
    theme_void()
  plot=ggarrange(plot1,plot2,plot3,ncol=1,common.legend = T,legend = "bottom")
  plot = annotate_figure(plot, top = text_grob(paste("Model: ", Hemisphere, sep=""), color = "black", face = "bold", size = 14)) # give figs a title
  return(plot)
}
dfp=function(data,Hemisphere){
  data$hemi=ifelse(grepl("lh_",data$Parameter)==T,"left","right")
  data = data[order(data$Parameter),]
  vols=data[grepl("_volume",data$Parameter),]
  thick=data[grepl("_thickness",data$Parameter),]
  surf=data[grepl("_area",data$Parameter),]
  vols$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(gsub("_volume","",vols$Parameter),3)][1:34]
  thick$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(gsub("_thickness","",thick$Parameter),3)][1:34]
  surf$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(gsub("_area","",surf$Parameter),3)][1:34]
  #
  plot1=ggplot(vols) + geom_brain(atlas = dk,aes(fill = df_std),color="black")+
    scale_fill_viridis_c(option = "cividis", direction = -1)+
    #scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Feature Importance of Volumes", sep="")) + labs(fill='df') +
    theme_void()
  plot2=ggplot(thick) + geom_brain(atlas = dk,aes(fill = df_std),color="black")+
    scale_fill_viridis_c(option = "cividis", direction = -1)+
    #scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Feature Importance of Cortical Thickness", sep="")) + labs(fill='df') +
    theme_void()
  plot3=ggplot(surf) + geom_brain(atlas = dk,aes(fill = df_std),color="black")+
    scale_fill_viridis_c(option = "cividis", direction = -1)+
    #scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Feature Importance of Surface Area", sep="")) + labs(fill='df') +
    theme_void()
  plot=ggarrange(plot1,plot2,plot3,ncol=1,common.legend = T,legend = "bottom")
  plot = annotate_figure(plot, top = text_grob(paste("Model: ", Hemisphere, sep=""), color = "black", face = "bold", size = 14)) # give figs a title
  return(plot)
}
#
## for knockout feature importance
anatp=function(data,age_span){
  # First, edit df for plotting
  data$hemi=ifelse(grepl("lh_",data$Area)==T,"left","right")
  data = data[order(data$Area),]
  data$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(data$Area,3)][1:34]
  # Second do the plotting with the right settings
  p1=ggplot(data%>%select(-Mean_Difference,-X)) + geom_brain(atlas = dk,aes(fill = Cohens_d),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Age-span: ", age_span, sep="")) + labs(fill="Cohen's d") +
    theme_void()
  return(p1)
}
anatp2=function(data,age_span){
  data$hemi=ifelse(grepl("lh_",data$Area)==T,"left","right")
  data = data[order(data$Area),]
  data$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(data$Area,3)][1:34]
  #
  p2=ggplot(data %>% filter(Mean_Difference > -5000)) + geom_brain(atlas = dk,aes(fill = Mean_Difference),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(title=paste("Age-span: ", age_span, sep="")) + labs(fill='Years') +
    theme_void()
  return(p2)
}
anatpn=function(data){
  # First, edit df for plotting
  data$hemi=ifelse(grepl("lh_",data$Area)==T,"left","right")
  data = data[order(data$Area),]
  data$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(data$Area,3)][1:34]
  # Second do the plotting with the right settings
  p1=ggplot(data%>%select(-Mean_Difference,-X)) + geom_brain(atlas = dk,aes(fill = Cohens_d),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(fill="Cohen's d") +
    theme_void()
  return(p1)
}
anatp2n=function(data){
  data$hemi=ifelse(grepl("lh_",data$Area)==T,"left","right")
  data = data[order(data$Area),]
  data$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(data$Area,3)][1:34]
  #
  p2=ggplot(data %>% filter(Mean_Difference > -5000)) + geom_brain(atlas = dk,aes(fill = Mean_Difference),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(fill='Years') +
    theme_void()
  return(p2)
}
# Plot #####
# feature importance plotting ####
fdf[[2]]$Parameter=gsub("rh_","lh_",fdf[[2]]$Parameter)#correct left hemi labels (if wrongly labeled as right)
hem=c("Both", "Left", "Right")
fp2=fp=list()
for(i in 1:length(fdf)){
  fp[[i]]=impp(fdf[[i]],hem[i])
  fp2[[i]]=dfp(fdf[[i]],hem[i])
}
#set_scale_union(fp[[1]],fp[[2]],fp[[3]],scale=scale_fill_gradient2(low = "blue",mid = "white",high="red"))
p0=ggarrange(plotlist = fp, ncol=3, common.legend = T, legend = "bottom")
p00=ggarrange(plotlist = fp2, ncol=3, common.legend = T, legend = "bottom")
#
ggsave("results/feature_importance.pdf",p0,width=15,height = 5) # save feature importance plots
ggsave("results/feature_importance.jpg",p0,width=15,height = 5)
ggsave("results/splineshape.pdf",p00,width=15,height = 5) # save plots indicating spline shape
ggsave("results/splineshape.jpg",p00,width=15,height = 5)
c(mean(fdf[[1]]$df_std),mean(fdf[[2]]$df_std),mean(fdf[[3]]$df_std)) # check the mean degrees of freedom

#
# knockout plots ####
# across predictions
ap=cp=list() # age plot, cohen's d plot lists
for(i in 1:length(explain)){
  cp[[i]]=anatpn(explain[[i]])
  ap[[i]]=anatp2n(explain[[i]])
}
a=ggarrange(plotlist = ap,ncol=1,common.legend = T,legend = "bottom")
b=ggarrange(plotlist = cp,ncol=1,common.legend = T,legend = "bottom")
pknock=ggarrange(a,b)
ggsave("results/knock.pdf",pknock,width=7,height = 4)
ggsave("results/knock.jpg",pknock,width=7,height = 4)
#
#
# age stratefied
age = c("0-10", "10-20","20-30","30-40","40-50","50-60","60-70","70-80","80+") # age span vector for plotting
ap=cp=list() # age plot, cohen's d plot lists
# BOTH BRAIN AGES
for(i in 1:length(ldf)){
  cp[[i]]=anatp(ldf[[i]],age[i])
  ap[[i]]=anatp2(ldf[[i]],age[i])
}
p1=ggarrange(plotlist = cp, ncol=1, common.legend = T, legend = "bottom") # create figs
p2=ggarrange(plotlist = ap, ncol=1, common.legend = T, legend = "bottom")
#p1=annotate_figure(p1, top = text_grob("Knock-out importance in Cohen's d", color = "black", face = "bold", size = 10)) # give figs a title
#p2=annotate_figure(p2, top = text_grob("Knock-out importance in Years", color = "black", face = "bold", size = 10)) # give figs a title
#
p0.1 = ggarrange(p1,p2,ncol=2)
p0.1 = annotate_figure(p0.1, top = text_grob("Both hemisphere brain age", color = "black", face = "bold", size = 10)) # give figs a title
ggsave("results/knockout_imp.pdf",p0.1,width=7,height = 10)
ggsave("results/knockout_imp.jpg",p0.1,width=7,height = 10)
#
# RIGHT BRAIN AGE ONLY
for(i in 1:length(ldfR)){
  cp[[i]]=anatp(ldfR[[i]],age[i])
  ap[[i]]=anatp2(ldfR[[i]],age[i])
}
p1R=ggarrange(plotlist = cp, ncol=1, common.legend = T, legend = "bottom") # create figs
p2R=ggarrange(plotlist = ap, ncol=1, common.legend = T, legend = "bottom")
p0.2 = ggarrange(p1R,p2R,ncol=2)
p0.2 = annotate_figure(p0.2, top = text_grob("Right hemispheric brain age", color = "black", face = "bold", size = 10)) # give figs a title
ggsave("results/knockout_impR.pdf",p0.2,width=7,height = 10)
ggsave("results/knockout_impR.jpg",p0.2,width=7,height = 10)
#
#
# LEFT BRAIN AGE ONLY
for(i in 1:length(ldfL)){
  cp[[i]]=anatp(ldfL[[i]],age[i])
  ap[[i]]=anatp2(ldfL[[i]],age[i])
}
p1L=ggarrange(plotlist = cp, ncol=1, common.legend = T, legend = "bottom") # create figs
p2L=ggarrange(plotlist = ap, ncol=1, common.legend = T, legend = "bottom")
#p1L=annotate_figure(p1L, top = text_grob("Knock-out importance in Cohen's d", color = "black", face = "bold", size = 14)) # give figs a title
#p2L=annotate_figure(p2L, top = text_grob("Knock-out importance in Years", color = "black", face = "bold", size = 14)) # give figs a title
p1L=ggarrange(plotlist = cp, ncol=1, common.legend = T, legend = "bottom") # create figs
p2L=ggarrange(plotlist = ap, ncol=1, common.legend = T, legend = "bottom")
p0.3 = ggarrange(p1L,p2L,ncol=2)
p0.3 = annotate_figure(p0.3, top = text_grob("Left hemispheric brain age", color = "black", face = "bold", size = 10)) # give figs a title
# ggsave("results/knockout_impL.pdf",p0.3,width=7,height = 10)
# ggsave("results/knockout_impL.jpg",p0.3,width=7,height = 10)
#
#
# Put them all together
p0.4=ggarrange(p0.1,p0.2,p0.3,ncol=3)
p0.4 = annotate_figure(p0.4, top = text_grob("Knock-out importance", color = "black", face = "bold", size = 14))
#
ggsave("results/stratified_knockout.pdf",p0.4,width=20,height = 10)
ggsave("results/stratified_knockout.jpg",p0.4,width=20,height = 10)
#
#
#
#
#
#
#
print("ALL DONE.")