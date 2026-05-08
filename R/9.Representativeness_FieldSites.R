#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Exampining representativeness of field sites relative to area where model predictions were made
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit May 8, 2026
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

rm(list=ls())
library(ggplot2)
library(terra)
library(ggeffects)
library(glmmTMB)
library(cowplot)
library(performance)
options(scipen = 999)
library(boot)

#read in field data
seed.df = read.csv(#Your path here
  "/Data/Seedling_Data_Processed_ClimDist_March2326.csv")
#Same data processing as in all scripts.
seed.df = seed.df[seed.df$Total_seed_ha<400000,]
seed.df[seed.df$Interval=="Long","Distance_Reburnedge_m"]=0
seed.df$Interval=as.factor(seed.df$Interval)
seed.df$Broad_seed_ha = seed.df$Total_seed_ha - seed.df$Con_seed_ha
seed.df$Pct_Con = seed.df$Con_seed_ha/seed.df$Total_seed_ha*100
seed.df[seed.df$Total_seed_ha==0,"Pct_Con"]=0
seed.df$ConRegen = NA
seed.df[seed.df$Pct_Con>=61,"ConRegen"]=1 #Conifers dominate post-fire
seed.df[is.na(seed.df$ConRegen),"ConRegen"]=0 #Conifers do not dominate post-fire

#Read in modelled trajectories
regen.df = read.csv("/Data/Regen_Trajectories_Data.csv")

seed.df$Obs_Preds = "Obs"
regen.df$Obs_Preds = "Preds"

seed.df=seed.df[,c("CMI","GDD","Distance_Reburnedge_m","Obs_Preds")]
regen.df = regen.df[,c("CMI","GDD","Distance_Reburnedge_m","Obs_Preds")]

#What is the range of CMI at the field sites
min(seed.df$CMI);max(seed.df$CMI)#-6.11 19.07
#What is the range of CMI at the pixels where predictions were made?
cmi.dens = ggplot(regen.df,aes(x=CMI))+
  geom_histogram(aes(y = after_stat(count / sum(count))),fill="pink",colour="black",bins=35) +
  geom_vline(xintercept=-6.11,linewidth=1)+
  geom_vline(xintercept=19.07,linewidth=1)+
  scale_y_continuous(limits=c(0,0.35),expand=c(0,0))+
  scale_x_continuous(breaks=seq(-15,90,by=10),expand=c(0,0))+
  labs(y = "Frequency")+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
cmi.dens

range(seed.df$CMI)#-6.11, 19.07
quantile(regen.df$CMI,0.905)#19.45
quantile(regen.df$CMI,0.025)#-6.31
nrow(regen.df[regen.df$CMI>=-6.11&regen.df$CMI<=19.07,])/nrow(regen.df)*100#87.7% of CMI values
#In the National data fall within the range of CMI values in the field data.
cmi.unique = unique(regen.df$CMI)
length(cmi.unique[cmi.unique>=-6.11&cmi.unique<=19.07])/length(cmi.unique)*100#Weirdly, also 87.7%

#Same for GDD 
#Range of GDD in the field data
range(seed.df$GDD)#612, 1154
#Vs. Range of GDD at pixels
gdd.dens = ggplot(regen.df,aes(x=GDD))+
  geom_histogram(aes(y = after_stat(count / sum(count))),fill="pink",colour="black",bins=30) +
  geom_vline(xintercept=612,linewidth=1)+
  geom_vline(xintercept=1154,linewidth=1)+
  scale_y_continuous(limits=c(0,0.21),expand=c(0,0))+
  scale_x_continuous(breaks=seq(380,1605,by=100),expand=c(0,0))+
  labs(y = "Frequency")+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
gdd.dens

range(seed.df$GDD)#612, 1154
quantile(regen.df$GDD,0.91)#1154
quantile(regen.df$GDD,0.01)#689.7
nrow(regen.df[regen.df$GDD>=612&regen.df$GDD<=1154,])/nrow(regen.df)*100 #90.6% of GDD values
#In the National data fall within the range of GDD values in the field data.
gdd.unique = unique(regen.df$GDD)
length(gdd.unique[gdd.unique>=612&gdd.unique<=1154])/length(gdd.unique)*100 #89.8% of values

#Same thing for distance to reburn edge
range(seed.df$Distance_Reburnedge_m)#0, 362.5
dist.dens = ggplot(regen.df,aes(x=Distance_Reburnedge_m))+
  geom_histogram(aes(y = after_stat(count / sum(count))),fill="pink",colour="black",bins=50) +
  geom_vline(xintercept=362.5,linewidth=1)+
  #geom_vline(xintercept=250,linewidth=1)+
  scale_y_continuous(limits=c(0,0.4),expand=c(0,0))+
  scale_x_continuous(breaks=seq(0,3940,by=250),expand=c(0,0))+
  labs(y = "Frequency",x="Distance to reburn edge (m)")+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dist.dens

range(seed.df$Distance_Reburnedge_m)#0, 362.5
quantile(regen.df$Distance_Reburnedge_m,0.82)#360.6
quantile(regen.df$Distance_Reburnedge_m,0)#100
nrow(regen.df[regen.df$Distance_Reburnedge_m<=362.5,])/nrow(regen.df)*100#83.5% of Distance values
#In the National data fall within the range of Distance values in the field data.
dist.unique = unique(regen.df$Distance_Reburnedge_m)
(length(dist.unique[dist.unique<362.5]))/(length(dist.unique))*100#1.7% of values. ERk. But this has a crazy tail.

densplots = plot_grid(cmi.dens,gdd.dens,dist.dens,labels=c("a","b","c"),nrow=3)
densplots

png(filename = #Your path here
      "/Figures/Extended Data/RealDat_Range_Represent_March92026.png",res=300,
    units = "cm", width = 15, height = 20)
densplots
dev.off()

#Alright, let's take another look at the stacking of these
#For example, are any pixels in locations where all three predictors don't appear in the field data?

# setwd(#Your path here
#   "/Data/Mod_Pred_Rasts/")
# prast = list.files(pattern="Predict")
# prast <- sprc(prast)
# prast = mosaic(prast)
# cmi=prast[[1]]
# cmi[cmi$CMI>19.07]=9999
# cmi[cmi$CMI<(-6.11)]=9999
# cmi[cmi!=9999]=0
# cmi[cmi$CMI==9999]=1
# plot(cmi)
# 
# gdd=prast[[2]]
# gdd[gdd$GDD>1154]=9999
# gdd[gdd$GDD<612]=9999
# gdd[gdd!=9999]=0
# gdd[gdd$GDD==9999]=1
# plot(gdd)
# 
# dist = prast[[3]]
# dist[dist$Distance_Reburnedge_m>362.5]=9999
# dist[dist!=9999]=0
# dist[dist$Distance_Reburnedge_m==9999]=1
# 
# clim = cmi+gdd
# conf = cmi+dist+gdd
# c2 = conf-clim
# freq(clim)
# freq(c2)
# freq(conf)#75.5% of areas are 100% within the range of the data
# #23.4% has 1 variable outside of the observed range #Of this 20%, 70% is from distance to edge, not from climate...
# #1.1% has 2 variables outside of the observed range
# #If you want to, you can plot this, or save it as a map of model prediction confidence
# plot(conf)
# writeRaster(conf,
#             #Your path here
#             "/Data/Mod_Pred_Rasts/Confidence_OutofRange.tif")
