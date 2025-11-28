#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Fitting models to post-fire tree regeneration field data to predict stem density and % conifer
#Test data of predicted stem density and conifer dominance is available at:
#https://github.com/EllWhitman/Canada_2023_Reburn/tree/main/data
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Nov 25, 2025
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

rm(list=ls())
library(terra)
library(ggplot2)
library(glmmTMB)
library(plyr)
library(cutpointr)
options(scipen = 999)

#Read in predicted post-reburn stem density
pregen = rast(#Your path here
  "/Data/Regen_StemDens_Predicted.tif")
#Read in grid of predicted post-fire conifer dominance
pcon = rast(#Your path here
  "/Data/Conifer_Dominance_Likelihood_Predicted.tif")

#Convert to dataframes (each pixel is a row), combine
pcon.df = as.data.frame(pcon)
pregen.df = as.data.frame(pregen)
cmi.df = as.data.frame(clim[["CMI"]])#My computer doesn't have enough memory to do this all as one at the national scale, other peoples' probably do.
gdd.df = as.data.frame(clim[["GDD"]])
dist.df = as.data.frame(clim[["Distance_Reburnedge_m"]])
rm(clim)

#Turn into one dataframe and rename columns
regen.df = cbind(pcon.df,pregen.df$X1)
regen.df = cbind(regen.df,dist.df$Distance_Reburnedge_m)
regen.df = cbind(regen.df,cmi.df$CMI)
regen.df = cbind(regen.df,gdd.df$GDD)

#Rename and save resulting dataframe of all pixels
names(regen.df) = c("PConifer","Total_Stems_ha","Distance_Reburnedge_m","CMI","GDD")
write.csv(regen.df,#Your Path here
          "/Data/Regen_Trajectories_Data.csv",row.names=T)

rm(list=ls())
gc()
regen.df = read.csv("/Users/ewhitman/Documents/Projects/Reburn_Occurence_2023/Data/Regen_Trajectories_Data.csv")

#Using ruleset identify predicted post-fire recovery trajectory of each pixel
#1000 stems/ha = sparse
#<.56 likelihood of conifer dominance = broadleaf
regen.df$Trajectory = NA
regen.df[regen.df$Total_Stems_ha<=1000,"Trajectory"] = "Sparse"
regen.df[regen.df$Total_Stems_ha>1000&regen.df$PConifer<0.56,"Trajectory"]="Broadleaf"
regen.df[is.na(regen.df$Trajectory),"Trajectory"]="Conifer"
table(regen.df$Trajectory)

#Create class of sparse vs. not
regen.df$Sparse_Class = NA
regen.df[regen.df$Total_Stems_ha<=1000,"Sparse_Class"] = 1
regen.df[regen.df$Total_Stems_ha>1000,"Sparse_Class"] = 0

set.seed(123)
regen.samp = sample(row.names(regen.df),size=1053)#random .1% sample
regen.samp = regen.df[row.names(regen.df)%in%regen.samp,]

#Calculate optimal cutpoint of distance from reburn edge after which regen is sparse
point=cutpointr(regen.samp,Distance_Reburnedge_m,Sparse_Class,tol_metric = 0.05,method = maximize_metric, 
                metric = youden,break_ties = min)
t.thresh = point$optimal_cutpoint #Optimal cutpoint = 300
t.thresh
plot(point)
write.csv(regen.df,#Your path here
          "/Data/Regen_Trajectories_Data.csv",row.names=F)
rm(list=ls())
gc()

#Read the new trajectories file back in
regen.df = read.csv(#Your path here
  "/Data/Regen_Trajectories_Data.csv")

#Create ECDF plots of post-fire trajectories as a function of predictor variables
cmi.p = ggplot(regen.df,aes(CMI,colour=Trajectory))+
  stat_ecdf(geom="step",linewidth=0.8)+
  scale_colour_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide=F)+
  ylab("")+
  xlab("Climate Moisture Index (1991 - 2020)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill='transparent'),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.y= element_text(colour="black",
                                  size=11),,
        axis.text.x= element_text(colour="black",
                                  size=11))
cmi.p

gdd.p = ggplot(regen.df,aes(GDD,colour=Trajectory))+
  stat_ecdf(geom="step",linewidth = 0.8)+
  scale_colour_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
  xlab("Growing Degree Days (1991 - 2020)")+
  ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill='transparent'),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.y= element_text(colour="black",
                                  size=11),
        axis.text.x= element_text(colour="black",
                                  size=11))
gdd.p


dist.p = ggplot(regen.df,aes(Distance_Reburnedge_m,colour=Trajectory))+
  stat_ecdf(geom="step",linewidth=.8)+
  scale_colour_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
  #geom_vline(xintercept = 300)+
  xlab("Distance to reburn edge (m)")+
  ylab("")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill='transparent'),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.y= element_text(colour="black",
                                  size=11),
        axis.text.x= element_text(colour="black",
                                  size=11))
dist.p

library(cowplot)
plots.comb = cowplot::plot_grid(cmi.p,gdd.p,dist.p,
                                nrow = 2,
                                align = "v",
                                labels=c("a","b","c"))
plots.comb

png(filename = #Your path here
      "/Figures/ECDFS_Trajectory_Drivers.png",res=300,
    units = "cm", width = 21, height = 15)
plots.comb
dev.off()
rm(list=ks())
gc()
