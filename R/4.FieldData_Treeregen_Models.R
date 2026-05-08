#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Fitting models to post-fire tree regeneration field data to predict stem density and % conifer
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit May 8, 2026
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

rm(list=ls())
library(ggplot2)
library(terra)
library(ggeffects)
library(glmmTMB)
library(cutpointr)
library(pscl)
library(geotargets)
library(cowplot)
library(performance)
library(ROSE)
options(scipen = 999)

#Data availability note: not all tree recruitment data used for these models is publicly available, therefore this file is not provided
#Similar models can be fitted to consolidated field datasets that are publicly available e.g: https://doi.org/10.3334/ORNLDAAC/1955
#The research team identified and downloaded polygons of fire perimeters for all overlapping fires examined in these studies. 
#Fire perimeter polygons were sourced from https://cwfis.cfs.nrcan.gc.ca/ (Canada) and https://services3.arcgis.com/T4QMspbfLg3qTGWY/ArcGIS/rest/services/Alaska_Fire_History_Up_To_2024/FeatureServer (Alaska)
#Climate data was sourced from https://ftp.maps.canada.ca/pub/nrcan_rncan/Climate-archives_Archives-climatologiques/NAM_30_Year_Averages/1991-2020/
#Field site status was allowed to override mapping errors, assuming that field visits were more reliable than GIS
#E.G., if a field site was called long interval or unburned, but it was slightly inside of recent a fire perimeter, this was assumed to be a mapping error.
#Two sites that were far from any detectable burned area were excluded.
#The output of these models is available for download at:#TBD
#A spatial output 'test' dataset is found in the data folder of this repo.

seeds = vect(#Your path here. Shapefile of field measures of tree species and density at recently burned field sites
  "/Data/Boreal_SeedlingRegen_Reburns_AllsitesMarch2326.shp")
rebs = vect(#Your path here
  "/Data/Fielddata_Reburn_Perims.shp")#Reburned area polygons for FIELD SITES ONLY. Not the dataset of all reburns 1992-2023
seeds = project(seeds, crs(rebs))#Confirm same projection
gdd = rast(#Your path here.
  "/Data/mly60arcsecond_1991-2020/sg60_10.tif")#Growing degree days normals
cmi = rast(#Your path here.
  "/Data/mly60arcsecond_1991-2020/cmi60_sum.tif")#CMI normals

template = rast(ext(rebs),nlyrs=1,resolution = 30)#Make a template for rasterizing reburned area
crs(template) = crs(rebs)
rebs.rast = rasterize(rebs,template,background=2)
rebs.rast[rebs.rast!=2]=NA

rebdist = distance(rebs.rast,unit='m')#Distance to reburn edge for each pixel, in metres
rebdist = mask(rebdist,rebs)#Mask to reburned areas only
distseed = extract(rebdist,sites)#Assign a distance to edge for each field site
sites$Distance_Reburnedge_m = distseed$layer

#Export the distance raster in case you need it again
writeRaster(rebdist,
            #Your path here.
            "/Data/Reburn_Edge_Dist_m.tif")

#Add other climate predictors (cmi and gdd)
sites = project(sites,crs(cmi))
cmisites = extract(cmi,sites)
gddsites = extract(gdd,sites)
sites$CMI = cmisites$cmi60_sum
sites$GDD = gddsites$sg60_10

seed.df = as.data.frame(sites)#Turn the whole thing into a data frame for analysis

seed.df$Interval = as.factor(seed.df$Interval)
seed.df$Total_seed_ha = round(seed.df$total_sdln*10000,digits=0)#Convert stems/m2 into stems/ha
seed.df$Con_seed_ha = round(seed.df$conifer_sd*10000,digits=0)
#Save your work
write.csv(seed.df,
          #Your path here
          "/Data/Seedling_Data_Processed_ClimDist_March2326.csv")
rm(list=ls())#Let's clean this mess up.

#Okay. Let's get on with the modelling
seed.df = read.csv(#Your path here
  "/Data/Seedling_Data_Processed_ClimDist_March2326.csv")
seed.df = seed.df[seed.df$Total_seed_ha<400000,]#Removing outlier site with too much leverage (99th percentile seedling count)
#These sites tended to make the models love to way overpredict.
seed.df[seed.df$Interval=="Long","Distance_Reburnedge_m"]=0#If long interval was measured in the field, distance to edge is 0

seed.df$Interval=as.factor(seed.df$Interval)
seed.df$Broad_seed_ha = seed.df$Total_seed_ha - seed.df$Con_seed_ha
seed.df$Pct_Con = seed.df$Con_seed_ha/seed.df$Total_seed_ha*100
seed.df[seed.df$Total_seed_ha==0,"Pct_Con"]=0
seed.df$ConRegen = NA
seed.df[seed.df$Pct_Con>=61,"ConRegen"]=1 #Conifers dominate post-fire
seed.df[is.na(seed.df$ConRegen),"ConRegen"]=0 #Conifers do not dominate post-fire

#Short interval sites have significantly lower seedling densities than longer interval
median(seed.df[seed.df$Interval=="Long","Total_seed_ha"]) #15500
median(seed.df[seed.df$Interval=="Short","Total_seed_ha"]) #4553.5
median(seed.df[seed.df$Interval=="Long","Con_seed_ha"]) #14575
median(seed.df[seed.df$Interval=="Short","Con_seed_ha"]) #300
#Are the differences significant?
wilcox.test(seed.df[seed.df$Interval=="Long","Total_seed_ha"],seed.df[seed.df$Interval=="Short","Total_seed_ha"])
wilcox.test(seed.df[seed.df$Interval=="Long","Con_seed_ha"],seed.df[seed.df$Interval=="Short","Con_seed_ha"]) #This is mainly explained by fewer conifers
wilcox.test(seed.df[seed.df$Interval=="Long","Broad_seed_ha"],seed.df[seed.df$Interval=="Short","Broad_seed_ha"]) #But there's fewer broadleaves as well.


#How many tree seedlings can establish after fire, as a function of GDD+CMI+Distance to seed source?
#Comparing model forms with nested structure
#Non-significant variables (a = 0.05) were removed from consideration
pseed = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
               ziformula=~GDD+Distance_Reburnedge_m,
               family=poisson(),
               data=seed.df)
pseed2 = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
             data=seed.df,
             family=poisson())
pseed3 = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post,
                 data=seed.df,
                 family=poisson())
check_overdispersion(pseed)#The models that are not zero-inflated are overdispersed.
check_overdispersion(pseed2)
check_overdispersion(pseed3)

AIC(logLik(pseed))
AIC(logLik(pseed2))
AIC(logLik(pseed3))#AIC of Pseed with random effect and ZI is lowest!

pR2(pseed)#McFadden's pR2 is similar for the ZI model and the poisson model. Non-mixed effect model is poor
pR2(pseed2)
pR2(pseed3)

rm(pseed2,pseed3)
#Visualizing the tree recruitment density model
seed.effplot = predict_response(pseed, terms = c("Distance_Reburnedge_m[all]","CMI[-5,10,20]","GDD[600,800,1000]","years_post[10]"),
                        type = "zero_inflated",
                        ci_level=0.95)

seed.effplot=plot(seed.effplot)+
  scale_colour_manual(values=c('#FFAAA1','#C2ADC5','#7BADD1'))+
  scale_fill_manual(values=c('#FFAAA1','#C2ADC5','#7BADD1'))+
  geom_hline(yintercept=1000)+
  geom_hline(yintercept=5000,linetype="dashed")+
  labs(y = "Tree seedling stems/ha",x = "Distance to reburn edge (m)",title="")+
  scale_y_continuous(breaks=c(0,20000,40000,60000,80000,100000,120000,140000),
                   label = c("0","20,000","40,000","60,000","80,000","100,000","120,000","140,000"))+
  theme(legend.position="none")

#Will those seedlings be conifers?
#Creating a more balanced sample (the raw data has 74% of sites dominated by conifers post-fire)
#Using a combination of random over and undersampling
set.seed(1234)
seed.over = ovun.sample(ConRegen~., data = seed.df, method = "both", N = 200)
table(seed.over$data$ConRegen)#Now there are 97 broadleaf, and 103 conifer-dominated sites

#Model of likelihood of post-fire conifer dominance
#GDD is not significant, possibly because conifers are quite cold tolerant?
pcon = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m,
             data=seed.over$data,
             family=binomial)
pcon2 = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m+(1|project_id),
               data=seed.over$data,
               family=binomial)

AIC(logLik(pcon))
AIC(logLik(pcon2))#AIC is improved without the random effect. User the simpler model.
check_overdispersion(pcon)
check_overdispersion(pcon2)

pR2(pcon)#McFadden's pR2 is identical between the two models
pR2(pcon2)
rm(pcon2)

#Visualize the conifer recruitment likelihood model
pcon.plot <- predict_response(pcon, terms = c("Distance_Reburnedge_m[all]","CMI[-5,10,20]"),
                        type = "fixed",
                        ci_level=0.95)

pcon.plot=plot(pcon.plot)+
  scale_colour_manual(values=c('#FFAAA1','#C2ADC5','#7BADD1'))+
  scale_fill_manual(values=c('#FFAAA1','#C2ADC5','#7BADD1'))+
  labs(y = "Conifer dominance likelihood",x = "Distance to reburn edge (m)",title="")
  


plots.comb = cowplot::plot_grid(seed.effplot,pcon.plot,
                                labels=c("a","b"))
plots.comb

#Save the figure
png(filename = #Your Path Here
      "/Figures/Extended Data/ED4_MarginalEffects_Seedling_Conifer_Model.png",res=300,
    units = "cm", width = 20, height = 10)
plots.comb
dev.off()
#Cleaning up as we shift to prediction
rm(pcon.plot,plots.comb,seed.effplot)
gc()

#Calculate a breakpoint to determine what predicted probability is best suited to define conifer dominance
#At field sites
pred = predict(pcon,type="response",newdata = seed.df,re.form=NA)
#Probably should set seed, but I get the same answer every time anyway...
point=cutpointr(seed.df,pred,ConRegen,tol_metric = 0.05,method = maximize_metric, 
                metric = youden,break_ties = min,boot_runs=1000,use.midpoints=T)
point

t.thresh = point$optimal_cutpoint #Optimal cutpoint = 0.55
t.thresh

#Let's create spatial predictions from the models
#Import CMI grid and Reburn vector, remove small reburns, reproject all to EA Projection
cmi.20 = rast(#Your path here.
  "/Data/mly60arcsecond_1991-2020/cmi60_sum.tif")
t.crs = "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Lambert Azimuthal Equal Area\",\n        METHOD[\"Lambert Azimuthal Equal Area\",\n            ID[\"EPSG\",9820]],\n        PARAMETER[\"Latitude of natural origin\",45,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]]]"
rebs = vect(#Your path here.
  "/Data/Reburns_TimeSeries_Complete_19802024_20yr.shp")
cmi.20 = project(cmi.20,crs(t.crs),method="bilinear")
rebs = project(rebs,crs(t.crs));rm(t.crs)
rebs = rebs[rebs$YEAR==2023,]#2023 only
rebs$AREA_Ha = expanse(rebs,"ha")
rebs = rebs[rebs$AREA_Ha>=4,]
sum(rebs$AREA_Ha)#1053249 Alright we still have the right amount. That's good.

#Make a template for rasterizing climate data with a 1ha pixel size
template = rast(ext(rebs),nlyrs=1,resolution = 100)
crs(template) = crs(rebs)
rebs.rast = rasterize(rebs,template)

#Assign CMI (climate moisture index) values from 1991-2020 normals to each area reburned in 2023
cmi.reb.20 = resample(cmi.20,rebs.rast,method="bilinear")
cmi.reb.20 = terra::mask(cmi.reb.20,rebs.rast)#Reburns only
CMI = cmi.reb.20;rm(cmi.reb.20,cmi.20)

#Read in ggd as well
gdd.20 = rast(#Your path here.
  "/Data/mly60arcsecond_1991-2020/sg60_10.tif")
gdd.20 = project(gdd.20,crs(CMI),method="bilinear")
#Assign GDD (growing degree day) values from 1991-2020 normals to each area reburned in 2023
gdd.reb.20 = resample(gdd.20,rebs.rast,method="bilinear")
gdd.reb.20 = terra::mask(gdd.reb.20,rebs.rast)#Only reburns
GDD = gdd.reb.20;rm(gdd.reb.20,gdd.20)

#Read in large (>4ha) reburns with Canvec 250k waterbodies added 
#Waterbodies must be filled, otherwise you are telling the model predictions that water is a seed source
#As the distance to edge (of water) in the interior of a reburn is 0...
#I filled in the waterbodies in QGIS
rebs.rast2 = vect(#Your Path here
  "/Data/Reburns_2023_Waterbodies_250k.shp")
rebs.rast2 = rasterize(rebs.rast2,template,background=2)
rebs.rast2[rebs.rast2!=2]=NA
rebdist = distance(rebs.rast2,unit='m')#Re-calculate distance to reburn edge, this time for the full 2023 dataset to use as a predictor
rebdist = terra::mask(rebdist,GDD)
gc()

preddat = c(CMI,GDD,rebdist)#We now have a raster stack of all predictors used by the models (except stand age but that will be held constant)
names(preddat) = c("CMI","GDD","Distance_Reburnedge_m")
rm(GDD,CMI,rebdist,rebs,rebs.rast,rebs.rast2,template)

#This dataset is quite a lot too big to be workable processing on my machine.
#To address this I tiled it, and processed tiles individually.
tilez = rast(ncols = 3,nrows=3,crs=crs(preddat),extent = ext(preddat))
#plot(tilez)
ptile = makeTiles(preddat,tilez,
                  #Your path here
                  "/Data/Mod_Pred_Rasts/Seed_Predict_Dat_.tif",
                  overwrite = T)
rm(tilez,preddat,ptile)
gc()
rm(point,seed.df,seed.over,pcon.plot,plots.comb,pred,ptile,seed.effplot,t.thresh,tilez)

pth = #Your path here. Wherever you want to save all your predictors and model outputs
  "/Data/Mod_Pred_Rasts/"
i = 1

#In order to predict in a loop the predict needs to be its own function
rfun <- function(mod, dat, ...) {
  library(glmmTMB)
  predict(mod, dat, ...)
}

#Save the predictor data in tiled form
#This runs for several hours, on my machine, so be aware
#The regen density are particularly slow, I think due to the range of values (0 to 20K+)
#Model outputs are predicted at a constant timestep of 10 years post-fire
#Model predictions are of fixed effects only
#Currently written to use 7 cores
for(i in 1:9){
pdat = rast(paste(pth,"Seed_Predict_Dat_",i,".tif",sep=""))
regen.pred = terra::predict(pdat,pseed,type="response",const = data.frame(years_post=10),
                           cores = 7, cpkgs = "glmmTMB",re.form=NA)
con.pred = terra::predict(pdat,pcon,type="response",
                        cpkgs = "glmmTMB",
                       cores = 7, re.form=NA)

writeRaster(regen.pred,paste(pth,"Seed_Regen_Dens_",i,".tif",sep=""),overwrite=T)
writeRaster(con.pred,paste(pth,"Conifer_Prob_",i,".tif",sep=""),overwrite=T)
print(paste("Done with tile ",i,sep=""))
print(Sys.time())
i=i+1
}

#Woof. That was a lot.
#Let's clean up, and then put things back together.
rm(list=ls())
gc()
setwd(#Your path here.
  "/Data/Mod_Pred_Rasts/")
regenlist = list.files(pattern="Regen")#All of the seedling density tiles
conlist = list.files(pattern="Conifer")#All of the conifer likelihood tiles
regenrc <- sprc(regenlist)
conrc <- sprc(conlist)
pregen = mosaic(regenrc)
pcon = mosaic(conrc)
#Save the mosaiced layers two rasters
writeRaster(pregen,
            #Your path here.
            "/Data/Regen_Prob_Predicted.tif",overwrite=T)
writeRaster(pcon,
            #Your path here.
            "/Data/Conifer_Dominance_Likelihood_Predicted.tif",overwrite=T)
rm(list=ls())#Thanks, computer, for your efforts. You deserve a rest.
