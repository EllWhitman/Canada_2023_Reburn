#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Shared only for review. See data availability note.
#Test data of model outputs is available for use from script 5 onwards
#https://github.com/EllWhitman/Canada_2023_Reburn/tree/main/Test_Data
#Fitting models to post-fire tree regeneration field data to predict stem density and % conifer
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Nov 25, 2025
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

Data availability note: not all tree recruitment data used for these models is publicly available, therefore this file is not provided
#Similar models can be fitted to consolidated field datasets that are publicly available e.g: https://doi.org/10.3334/ORNLDAAC/1955
#The research team identified and downloaded polygons of fire perimeters for all overlapping fires examined in these studies. 
#Fire perimeter polygons were sourced from https://cwfis.cfs.nrcan.gc.ca/ (Canada) and https://services3.arcgis.com/T4QMspbfLg3qTGWY/ArcGIS/rest/services/Alaska_Fire_History_Up_To_2024/FeatureServer (Alaska)
#Climate data was sourced from https://ftp.maps.canada.ca/pub/nrcan_rncan/Climate-archives_Archives-climatologiques/NAM_30_Year_Averages/1991-2020/
#Field site status was allowed to override mapping errors, assuming that field visits were more reliable than GIS
#E.G., if a field site was called long interval or unburned, but it was slightly inside of recent a fire perimeter, this was assumed to be a mapping error
#The national output of these models will be made available on the Canadian Open Government Data Portal.
#A subset for use with the other scripts in this repo is available at: https://github.com/EllWhitman/Canada_2023_Reburn/tree/main/Test_Data

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
#E.G., if a field site was called long interval or unburned, but it was slightly inside of recent a fire perimeter, this was assumed to be a mapping error
#The output of these models is available for download at:

seeds = vect(#Your Path here. Shapefile of field measures of tree species and density at recently burned field sites
  "/Data/Boreal_SeedlingRegen_Reburns_AllsitesJun25.shp")
rebs = vect(#Your Path here. Shapefile of reburn overlaps for fires of interest, created using the first two scripts.
  #Be careful, if sites reburn multiple times, or if a fire is an old fire relative to one year, and the more recent
  #reburning fire in another. Can cause issues with the distance to edge calculation. 
  #If so, calculate distance to edge by year, instead of for all reburns at once
  "Fielddata_Reburn_Perims.shp")
seeds = project(seeds, crs(rebs))
gdd = rast(#Your path here
  "/Data/mly60arcsecond_1991-2020/sg60_10.tif")
cmi = rast(#Your path here
  "/Data/mly60arcsecond_1991-2020/cmi60_sum.tif")

template = rast(ext(rebs),nlyrs=1,resolution = 30)#Make a template for rasterizing
crs(template) = crs(rebs)
rebs.rast = rasterize(rebs,template,background=2)
rebs.rast[rebs.rast!=2]=NA#Set reburns as NA, backround as 2
rebdist = distance(rebs.rast,unit='m')#Calculate distance to reburn edge
rebdist = mask(rebdist,rebs)#Mask to reburned areas only
distseed = extract(rebdist,sites)#Extract distance to reburn edge at field sites
sites$Distance_Reburnedge_m = distseed$layer

sites = project(sites,crs(cmi))
cmisites = extract(cmi,sites)#Extract CMI normals at field sites
gddsites = extract(gdd,sites)#Extract GDD normals at field sites
sites$CMI = cmisites$cmi60_sum
sites$GDD = gddsites$sg60_10

seed.df = as.data.frame(sites)#Conslidate it all into one data frame
head(seed.df)

seed.df$Interval = as.factor(seed.df$Interval)
seed.df$Total_seed_ha = round(seed.df$total_sdln*10000,digits=0)#Convert from stems/m2 to stems/ha
seed.df$Con_seed_ha = round(seed.df$conifer_sd*10000,digits=0)
write.csv(seed.df,#Your path here
          "/Data/Seedling_Data_Processed_ClimDist.csv")
rm(list=ls())

#Read in the data for model building
seed.df = read.csv(#Your path here
  "/Data/Seedling_Data_Processed_ClimDist.csv")
seed.df[seed.df$Interval=="Long","Distance_Reburnedge_m"]=0#Just in case there were any mapping errors (explained above)
#Set distance to reburn edge at long interval sites to 0m

seed.df$Interval=as.factor(seed.df$Interval)
seed.df$Broad_seed_ha = seed.df$Total_seed_ha-seed.df$Con_seed_ha
seed.df$Pct_Con = seed.df$Con_seed_ha/seed.df$Total_seed_ha*100 #Calculate % conifer in post-fire cohort
seed.df[seed.df$Total_seed_ha==0,"Pct_Con"]=0#Set % conifer to zero if no trees
seed.df$ConRegen = NA#Create binomial factor for conifer dominance
seed.df[seed.df$Pct_Con>=61,"ConRegen"]=1 #Conifers dominate post-fire
seed.df[is.na(seed.df$ConRegen),"ConRegen"]=0 #Conifers do not dominate post-fire

#Short interval sites have significantly lower seedling densities than longer interval
median(seed.df[seed.df$Interval=="Long","Total_seed_ha"]) #16125
median(seed.df[seed.df$Interval=="Short","Total_seed_ha"]) #5000
median(seed.df[seed.df$Interval=="Long","Con_seed_ha"]) #15425
median(seed.df[seed.df$Interval=="Short","Con_seed_ha"]) #400
wilcox.test(seed.df[seed.df$Interval=="Long","Total_seed_ha"],seed.df[seed.df$Interval=="Short","Total_seed_ha"])
wilcox.test(seed.df[seed.df$Interval=="Long","Con_seed_ha"],seed.df[seed.df$Interval=="Short","Con_seed_ha"]) #This is mainly explained by fewer conifers
wilcox.test(seed.df[seed.df$Interval=="Long","Broad_seed_ha"],seed.df[seed.df$Interval=="Short","Broad_seed_ha"]) #But there's fewer broadleaves as well.


#How many tree seedlings can establish after fire?
#Nested glmmms predicting total seedling density in long and short-interval sites with a random effect of project ID
pseed= glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
               ziformula=~GDD+Distance_Reburnedge_m,
               family=poisson(),
               data=seed.df)
pseed2 = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
             data=seed.df,
             family=poisson())
pseed3 = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post,
                 data=seed.df,
                 family=poisson())
AIC(logLik(pseed))
AIC(logLik(pseed2))
AIC(logLik(pseed3))#AIC of Pseed with random effect and ZI is lowest!

pR2(pseed)
pR2(pseed2)
pR2(pseed3)

rm(pseed2,pseed3)
#Produce effect plot of seedling density model
seed.effplot = predict_response(pseed, terms = c("Distance_Reburnedge_m[all]","CMI[-5,10,20]","GDD[600,800,1000]","years_post[10]"),
                        type = "fixed",
                        ci_level=NA)

seed.effplot=plot(seed.effplot)+
  scale_colour_manual(values=c('#FFAAA1','#C2ADC5','#7BADD1'))+
  geom_hline(yintercept=1000)+
  geom_hline(yintercept=5000,linetype="dashed")+
  labs(y = "Tree seedling stems/ha",x = "Distance to reburn edge (m)",title="")+
  scale_y_continuous(breaks=c(0,20000,40000,60000,80000),
                   label = c("0","20,000","30,000","40,000","80,000"))+
  theme(legend.position="none")

#Will those seedlings be conifers?
#GDD is not significant. Dropped from model
set.seed(123)
#Correct imbalance between conifer dominated vs. not with ROSE methodology
seed.over = ovun.sample(ConRegen~., data = seed.df, method = "both", N = 200)
table(seed.over$data$ConRegen)

#Model of likelihood of post-fire conifer dominance
pcon = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m,
             data=seed.over$data,
             family=binomial)
pcon2 = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m+(1|project_id),
               data=seed.over$data,
               family=binomial)

AIC(logLik(pcon))
AIC(logLik(pcon2))#AIC is improved without the random effect. User the simpler model.

pR2(pcon)
pR2(pcon2)
rm(pcon2)
#Produce effect plot of conifer dominance model
pcon.plot <- predict_response(pcon, terms = c("Distance_Reburnedge_m[all]","CMI[-5,10,20]"),
                        type = "fixed",
                        ci_level=NA)

pcon.plot=plot(pcon.plot)+
  scale_colour_manual(values=c('#FFAAA1','#C2ADC5','#7BADD1'))+
  labs(y = "Conifer dominance likelihood",x = "Distance to reburn edge (m)",title="")
  


plots.comb = cowplot::plot_grid(seed.effplot,pcon.plot,
                                labels=c("a","b"))
plots.comb

png(filename = #Your path here
      "/Figures/MarginalEffects_Seedling_Conifer_Model.png",res=300,
    units = "cm", width = 20, height = 10)
plots.comb
dev.off()
rm(pcon.plot,plots.comb,seed.effplot)

gc()
#Calculate optimal cutpoint to determine whether the model prediction should be considered likely conifer dominance
pred = predict(pcon,type="response",newdata = seed.df,re.form=NA)
point=cutpointr(seed.df,pred,ConRegen,tol_metric = 0.05,method = maximize_metric, metric = youden,break_ties = min)
t.thresh = point$optimal_cutpoint #Optimal cutpoint = 0.56
t.thresh


#Import CMI grid and Reburn vector, remove small reburns, reproject all to EA Projection
cmi.20 = rast(#Your path here
  "/Data/mly60arcsecond_1991-2020/cmi60_sum.tif")
t.crs = "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Lambert Azimuthal Equal Area\",\n        METHOD[\"Lambert Azimuthal Equal Area\",\n            ID[\"EPSG\",9820]],\n        PARAMETER[\"Latitude of natural origin\",45,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]]]"
rebs = vect(#Your path here. This file was created with the first two scripts
  "/Data/Reburns_TimeSeries_Complete_19822023_20yr.shp")
cmi.20 = project(cmi.20,crs(t.crs),method="bilinear")
rebs = project(rebs,crs(t.crs));rm(t.crs)
rebs = rebs[rebs$YEAR==2023,]#2023 only
rebs$AREA_Ha = expanse(rebs,"ha")
rebs = rebs[rebs$AREA_Ha>=4,]
sum(rebs$AREA_Ha)#1053249

#Make a template for rasterizing climate data with a 1ha pixel size
template = rast(ext(rebs),nlyrs=1,resolution = 100)
crs(template) = crs(rebs)
rebs.rast = rasterize(rebs,template)

#Assign CMI (climate moisture index) values from 1991-2020 normals to each area reburned in 2023
cmi.reb.20 = resample(cmi.20,rebs.rast,method="bilinear")
cmi.reb.20 = terra::mask(cmi.reb.20,rebs.rast)#Reburns only
CMI = cmi.reb.20;rm(cmi.reb.20,cmi.20)

gdd.20 = rast(#Your path here
  "/Data/mly60arcsecond_1991-2020/sg60_10.tif")
gdd.20 = project(gdd.20,crs(CMI),method="bilinear")
#Assign GDD (growing degree day) values from 1991-2020 normals to each area reburned in 2023
gdd.reb.20 = resample(gdd.20,rebs.rast,method="bilinear")
gdd.reb.20 = terra::mask(gdd.reb.20,rebs.rast)#Only reburns
GDD = gdd.reb.20;rm(gdd.reb.20,gdd.20)

#In GIS fill holes in reburns from waterbodies. Waterbodies that are alongside or entirely inside of a reburn
#do not have seeds for trees generated from them
#Read in large (>4ha) reburns with Canvec 250k waterbodies added 
rebs.rast2 = vect(#Your path here
  "/Reburns_2023_Waterbodies_250k.shp")
rebs.rast2 = rasterize(rebs.rast2,template,background=2)#Conver to raster
rebs.rast2[rebs.rast2!=2]=NA
rebdist = distance(rebs.rast2,unit='m')#Calculate distance to reburn edge in 2023 reburns
rebdist = terra::mask(rebdist,GDD)#Mask to reburned areas only
gc()

preddat = c(CMI,GDD,rebdist)#Combine the model predictor rasters in all 2023 reburns
names(preddat) = c("CMI","GDD","Distance_Reburnedge_m")
rm(GDD,CMI,rebdist,rebs,rebs.rast,rebs.rast2,template)#Clean up your workspace

tilez = rast(ncols = 3,nrows=3,crs=crs(preddat),extent = ext(preddat))#This data is otherwise way too big to export on my computer
#This will break Canada into 9 tiles, so that the code doesn't crash while exporting the rasters of predictor variables
#plot(tilez)
#Save tiles of predictor variable rasters in reburned areas
ptile = makeTiles(preddat,tilez,
                  #Your path here
                  "/Data/Mod_Pred_Rasts/Seed_Predict_Dat_.tif",
                  overwrite = T)
rm(tilez,preddat,ptile)
#clean up
gc()
rm(point,seed.df,seed.over,pcon.plot,plots.comb,pred,ptile,seed.effplot,t.thresh,tilez)

#Your path here. Where you saved the tiles
pth = "/Data/Mod_Pred_Rasts/"
i = 1

rfun <- function(mod, dat, ...) {
  library(glmmTMB)
  predict(mod, dat, ...)
}

#Predict stem density of post-fire recruitment at 10 years after a fire, and likelihood of conifer dominance
#In all 2023 reburned areas
for(i in 1:9){
pdat = rast(paste(pth,"Seed_Predict_Dat_",i,".tif",sep=""))
regen.pred = terra::predict(pdat,pseed,type="response",const = data.frame(years_post=10),
                            cores = 12, cpkgs = "glmmTMB",re.form=NA)
con.pred = terra::predict(pdat,pcon,type="response",
                          cpkgs = "glmmTMB",
                          cores = 12, re.form=NA)

writeRaster(regen.pred,paste(pth,"Seed_Regen_Dens_",i,".tif",sep=""),overwrite=T)
writeRaster(con.pred,paste(pth,"Conifer_Prob_",i,".tif",sep=""),overwrite=T)
print(paste("Done with tile ",i,sep=""))
i=i+1
}

rm(list=ls())
gc()

setwd(#Your path here
  "/Data/Mod_Pred_Rasts/")
regenlist = list.files(pattern="Regen")
conlist = list.files(pattern="Conifer")
regenrc <- sprc(regenlist)
conrc <- sprc(conlist)
pregen = mosaic(regenrc)#Mosaic the tiles back together and save the final predicted stem density and conifer likelihood rasters
pcon = mosaic(conrc)
writeRaster(pregen,"/Data/Regen_StemDens_Predicted.tif",overwrite=T)
writeRaster(pcon,"/Data/Conifer_Dominance_Likelihood_Predicted.tif",overwrite=T)
rm(list=ls())
