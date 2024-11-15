# Analyse multiverse outputs
#
# Max Korbmacher, 14 Nov 2024
#
# recap, number of models run:
choose(19,10)*2*2

# load stuff
library(dplyr)
library(sjstats)
#devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)
#
# Then load specific for the two predictions
#
#
#
# PASAT
df=read.csv("/Users/max/Documents/Local/MS/results/PASAT_effects.csv")
# get estimates
l=df %>% group_by(Names)%>%summarize(Beta_Md=median(Beta), Beta_MAD = mad(Beta), P_Md=median(P)) # estimate medians and MADs
a = props(df%>%group_by(Names),Beta>0) # estimate proportions of effects towards a single direction.
b = props(df%>%group_by(Names),Beta<0)
l = merge(l,a,by="Names") # merge data keeping the correct order
l = merge(l,b,by="Names")

ggforestplot::forestplot(
  df = l,
  name = Names,
  estimate = Beta_Md,
  se = Beta_MAD
)
#
#
#
#
# EDSS
df2=read.csv("/Users/max/Documents/Local/MS/results/EDSS_effects.csv")
# get estimates
edss=df2 %>% group_by(Names)%>%summarize(Beta_Md=median(Beta), Beta_MAD = mad(Beta), P_Md=median(P)) # estimate medians and MADs
a = props(df2%>%group_by(Names),Beta>0) # estimate proportions of effects towards a single direction.
b = props(df2%>%group_by(Names),Beta<0)
edss = merge(edss,a,by="Names") # merge data keeping the correct order
edss = merge(edss,b,by="Names")

ggforestplot::forestplot(
  df = edss,
  name = Names,
  estimate = Beta_Md,
  se = Beta_MAD
)
