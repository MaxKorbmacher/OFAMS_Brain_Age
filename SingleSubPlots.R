# Individual-level metric-by-region brain age dependencies
#
# Max Korbmacher 02. December 2024
#
# pkgs
library(dplyr)
# data
## both hemispheres
path = "/Users/max/Library/CloudStorage/OneDrive-HøgskulenpåVestlandet/Documents/Projects/Hemi2/Explain_model/single_sub_results"
filenames <- list.files(path, pattern="*.csv", full.names=TRUE)
ldf <- lapply(filenames, read.csv)
both = Reduce(full_join,ldf)
both # shows that correlations between numeric changes and brain age per metric and region are the same across the age span
#
explain=list(both)
# create metric-specific lists
vol = surf = thick = list()
for (i in 1:length(explain)){
  vol[[i]] = explain[[i]] %>% filter(grepl("volume",Region))
  vol[[i]]$Region = gsub("_volume","",vol[[i]]$Region)
  surf[[i]] = explain[[i]] %>% filter(grepl("area",Region))
  surf[[i]]$Region = gsub("_area","",surf[[i]]$Region)
  thick[[i]] = explain[[i]] %>% filter(grepl("thicknes",Region))
  thick[[i]]$Region = gsub("_thickness","",thick[[i]]$Region)
}
#
# Plot
anatpn=function(data){ # make plot function
  # First, edit df for plotting
  data$hemi=ifelse(grepl("lh_",data$Region)==T,"left","right")
  data = data[order(data$Region),]
  data$region = brain_regions(dk)[substring(brain_labels(dk),3) %in% substring(data$Region,3)][1:34]
  # Second do the plotting with the right settings
  p1=ggplot(data) + geom_brain(atlas = dk,aes(fill = Correlation),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    labs(fill="Spearman's rho") +
    theme_void()
  return(p1)
}
v=s=t=list() # age plot, cohen's d plot lists
for(i in 1:length(vol)){
  v[[i]]=anatpn(vol[[i]])
  s[[i]]=anatpn(surf[[i]])
  t[[i]]=anatpn(thick[[i]])
}
plot = ggarrange(plotlist = c(v,s,t), ncol = 1, nrow = 3, common.legend = T, legend = "bottom", 
                 labels = c("Volume","Area","Thickness"),
                 label.x = c(0.39,0.425,0.35),label.y = 1)
path = "/Users/max/Documents/Local/MS/results"
ggsave(paste(path,"/single_sub.pdf",sep=""),plot,width=6,height = 5)
ggsave(paste(path,"/single_sub.jpg",sep=""),plot,width=6,height = 5)
#
# Another way of plotting (yet running into trouble with legend)
# # create plots for right,left,both in one figure, per metric
# a=ggarrange(plotlist = v,ncol=1,common.legend = T,legend = "bottom")
# b=ggarrange(plotlist = s,ncol=1,common.legend = T,legend = "bottom")
# c=ggarrange(plotlist = t,ncol=1,common.legend = T,legend = "bottom")
# # label plots
# a = annotate_figure(a, top = text_grob("Volume", color = "black", face = "bold", size = 10))
# b = annotate_figure(b, top = text_grob("Surface Area", color = "black", face = "bold", size = 10))
# c = annotate_figure(c, top = text_grob("Thickness", color = "black", face = "bold", size = 10))
# ggarrange(a,b,c, nrow=1, common.legend = T)
