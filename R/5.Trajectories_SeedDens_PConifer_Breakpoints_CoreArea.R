#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Fitting models to post-fire tree regeneration field data to predict stem density and % conifer
#Analysing modelled post-fire vegetation community trajectories
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit May 8, 2026
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file


rm(list=ls())
library(terra)
library(ggplot2)
library(glmmTMB)
library(plyr)
library(cutpointr)
library(ROSE)
options(scipen = 999)

#Read in field-measured seedling data, and reproduce models of stem density and % conifer fitted in script 5
seed.df = read.csv(#Read in seedling and predictor data
  "/Data/Seedling_Data_Processed_ClimDist_March2326.csv")
seed.df[seed.df$Interval=="Long","Distance_Reburnedge_m"]=0
seed.df = seed.df[seed.df$Total_seed_ha<400000,]#Remove outliers, as previously
seed.df$Interval=as.factor(seed.df$Interval)
seed.df$Broad_seed_ha = seed.df$Total_seed_ha - seed.df$Con_seed_ha
seed.df$Pct_Con = seed.df$Con_seed_ha/seed.df$Total_seed_ha*100
seed.df[seed.df$Total_seed_ha==0,"Pct_Con"]=0
seed.df$ConRegen = NA
seed.df[seed.df$Pct_Con>=61,"ConRegen"]=1 #Conifers dominate post-fire
seed.df[is.na(seed.df$ConRegen),"ConRegen"]=0 #Conifers do not dominate post-fire

#Replicate same balanced sample for conifer dominance model
set.seed(1234)
seed.over = ovun.sample(ConRegen~., data = seed.df, method = "both", N = 200)
table(seed.over$data$ConRegen)

#Model of likelihood of post-fire conifer dominance
pcon = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m,
               data=seed.over$data,
               family=binomial)
seed.df$Con_Pred = predict(pcon,type="response",newdata = seed.df)

#What breakpoint best predicts post-fire conifer dominance?
point=cutpointr(seed.df,Con_Pred,ConRegen,tol_metric = 0.05,method = maximize_metric, 
                metric = youden,break_ties = min,boot_runs=1000,use.midpoints=T)
t.thresh = point$optimal_cutpoint #Optimal cutpoint = 0.55
t.thresh
plot(point)
boot_ci(x=point,variable = optimal_cutpoint,alpha = 0.1)#0.4 0.7


pseed = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
                ziformula=~Distance_Reburnedge_m+GDD,
                family=poisson(),
                data=seed.df)

#Read in predicted post-reburn stem density
pregen = rast("/Data/Regen_StemDens_Predicted.tif")
rm(seed.df)

#Reading in the predictor rasters
setwd(#Your path here.
  "/Data/Mod_Pred_Rasts/")
climlist = list.files(pattern="Seed_Predict")
climrc <- sprc(climlist)
clim = mosaic(climrc)
gc()

#Read in grid of predicted post-fire conifer dominance
pcon = rast("/Data/Conifer_Dominance_Likelihood_Predicted.tif")
rm(climrc,point,seed.over,climlist,t.thresh)

#Convert to rasters to data frames (each pixel is a row), combine
pcon.df = as.data.frame(pcon)
pregen.df = as.data.frame(pregen)
cmi.df = as.data.frame(clim[["CMI"]])#My computer doesn't have enough memory to do this all as one, other peoples' probably do.
gdd.df = as.data.frame(clim[["GDD"]])
dist.df = as.data.frame(clim[["Distance_Reburnedge_m"]])
rm(clim)

regen.df = cbind(pcon.df,pregen.df$X1)
regen.df = cbind(regen.df,dist.df$Distance_Reburnedge_m)
regen.df = cbind(regen.df,cmi.df$CMI)
regen.df = cbind(regen.df,gdd.df$GDD)

#Rename and save resulting df of predictors and predictions all pixels
names(regen.df) = c("PConifer","Total_Stems_ha","Distance_Reburnedge_m","CMI","GDD")
write.csv(regen.df,
          #Your path here
          "/Data/Regen_Trajectories_Data_Mar272026.csv",row.names=T)

rm(list=ls())
gc()
#Starting over with just the data frame
regen.df = read.csv(#your path here
  "Data/Regen_Trajectories_Data_March272026.csv")

#Using ruleset identify trajectory of each pixel
#1000 stems/ha = sparse
#<.55 likelihood of conifer dominance = broadleaf
regen.df$Trajectory = NA
regen.df[regen.df$Total_Stems_ha<=1000,"Trajectory"] = "Sparse"
regen.df[regen.df$Total_Stems_ha>1000&regen.df$PConifer<0.55,"Trajectory"]="Broadleaf"
regen.df[is.na(regen.df$Trajectory),"Trajectory"]="Conifer"
table(regen.df$Trajectory)
nrow(regen.df[regen.df$Total_Stems_ha<=100,])#Number of pixels (or hectares) very sparse

#Create class of sparse vs. not for defining core area
regen.df$Sparse_Class = NA
regen.df[regen.df$Trajectory=="Sparse","Sparse_Class"] = 1
regen.df[is.na(regen.df$Sparse_Class),"Sparse_Class"] = 0

set.seed(1234)
regen.samp = sample(row.names(regen.df),size=9526)#random ~1% sample
regen.samp = regen.df[row.names(regen.df)%in%regen.samp,]

#Calculate optimal cutpoint of distance from reburn edge after which regen is sparse
point=cutpointr(regen.samp,Distance_Reburnedge_m,Sparse_Class,tol_metric = 0.05,method = maximize_metric, 
                metric = youden,break_ties = min,boot_runs=1000,use.midpoints=T)
t.thresh = point$optimal_cutpoint #Optimal cutpoint = 300
t.thresh
plot(point)
boot_ci(x=point,variable = optimal_cutpoint,alpha = 0.1)#I am confusion. How did I get multiple values
#Before and now I'm getting just one value?


write.csv(regen.df,
          #Your path here
          "/Data/Regen_Trajectories_Data_March272026.csv",row.names=F)
rm(list=ls())
regen.df = read.csv(#Your path here
  "/Data/Regen_Trajectories_Data_March272026.csv")
table(regen.df$Trajectory)
nrow(regen.df[regen.df$Total_Stems_ha<=1000,])#Area sparse
nrow(regen.df[regen.df$Total_Stems_ha<=1000&regen.df$Distance_Reburnedge_m>300,])#Area sparse that is 
#Greater than 300m from reburn edge
table(regen.df$Trajectory)
table(regen.df[regen.df$GDD<800|regen.df$CMI<=0,"Trajectory"])#What is going on in cold or dry places?
142953/210671*100#Amount of sparse in cold/dry
nrow(regen.df[regen.df$CMI<=0,])#497312 ha in dry clims
nrow(regen.df[regen.df$GDD<=800,])#46948 ha in cold clims

#Area of each trajectory in dry and cool climates
clim.obs = c(392423+7810,142953)
#Proportion of each trajectory in full dataset
clim.exp = c(0.6531471+0.1468314,0.2000215)

chisq.test(x = clim.obs,p = clim.exp,rescale.p = TRUE)#P<0.001,  sparse trajectories are more likely below a CMI of 0 or GDD of 800

range(regen.df[regen.df$Trajectory=="Conifer","CMI"])#Minimum CMI for conifer dominance is 9.6
temp=regen.df[regen.df$Total_Stems_ha<=100,]
nrow(temp[temp$GDD<800|temp$CMI<=0,])#70% of sparse is in cold or dry
nrow(regen.df[regen.df$Total_Stems_ha<=100,])
median(temp$Distance_Reburnedge_m)#median distance to edge of very sparse 640.3
rm(temp)

range(regen.df[regen.df$Distance_Reburnedge_m>=360.6,"Total_Stems_ha"])
range(regen.df[regen.df$Distance_Reburnedge_m>=300,"Total_Stems_ha"])
range(regen.df[regen.df$Total_Stems_ha<=1000,"Distance_Reburnedge_m"])#Sparse 200 to 4000 m to edge
nrow(regen.df[regen.df$Total_Stems_ha<=100,])
range(regen.df[regen.df$Trajectory=="Conifer","Distance_Reburnedge_m"])#No conifer regen predicted beyond 360.6

nrow(regen.df[regen.df$CMI<=0,])#Area dry 497312
nrow(regen.df[regen.df$GDD<=800,]) #Area cold 46948
table(regen.df$Trajectory)

#Let's visualize how the vegetation trajectories play out along the predictors
#Climate Moisture Index
cmi.p = ggplot(regen.df,aes(CMI,colour=Trajectory))+
  stat_ecdf(geom="step",linewidth=0.8)+
  scale_colour_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide=F)+
  ylab("")+
  xlab("Climate moisture index (1991 - 2020)")+
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
cmi.p#Deal with it GGPLOT, I don't care whether it's 'false' or 'none'

#Growing Degree Days
gdd.p = ggplot(regen.df,aes(GDD,colour=Trajectory))+
  stat_ecdf(geom="step",linewidth = 0.8)+
  scale_colour_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
  xlab("Growing degree days (1991 - 2020)")+
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

#Distance to Reburn Edge
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

#Save that plot, baby!
png(filename = #Your path here
      '/Figures/Extended Data/ED5_ECDFS_Trajectory_Drivers.png',res=300,
    units = "cm", width = 21, height = 15)
plots.comb
dev.off()
rm(plots.comb,cmi.p,gdd.p,dist.p)
