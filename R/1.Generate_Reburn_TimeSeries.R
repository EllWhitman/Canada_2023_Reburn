#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
library(terra)

#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Feb 5, 2024
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file
#The ecozones file can be downloaded here: https://sis.agr.gc.ca/cansis/nsdb/ecostrat/gis_data.html
#The NBAC fire perimeter maps can be downloaded here:https://cwfis.cfs.nrcan.gc.ca/datamart
#If the path above ceases to function, access the NBAC maps via the main CWFIS page: https://cwfis.cfs.nrcan.gc.ca/
#Under the 'CWFIS datamart' tab

#Read in ecozones, and subset to boreal
zns = vect(#Your Path here
  "ecozones.shp")
#Yes there is a typo in Boreal PLain, and the capital L should remain.
#The analysis was originally produced with a slightly different ecozone map, used in carbon accounting applications
#where large ecozones (e.g. boreal shield) are split into east and west, but their outer extent should be the same
zns = zns[zns$ZONE_NAME%in%c("Taiga Shield","Taiga Plain","Boreal PLain","Boreal Shield",
                             "Taiga Cordillera","Boreal Cordillera","Hudson Plain"),]
#Read in fires, and clip to boreal
fires = vect(#Your path here
  "nbac_1972_2023_20240530.shp")
zns = project(zns,crs(fires)) #same projection
fires = fires[zns,]
writeVector(fires,
            #Your path here
            "Fires_Reburn_Merged_Clip.shp",overwrite=T)
rm(list=ls())
gc()


#Okay, we have our dataset for analysis
#Now let's identify some reburns!
fires = vect(#Your path here
  "Fires_Reburn_Merged_Clip.shp")
template = rast(ext(fires),nlyrs=1,resolution = 100)#Make a template for rasterizing
crs(template) = crs(fires)

#Run for all years, 1982 - 2023
yrs = seq(1982, 2023,by=1)
window = 20 #Start with a 20-yr window for simplicty, but 10-yr and 50-yr also of interest

for(i in 1:length(yrs)){
  print(paste("Running Year ", yrs[i],sep=""))
  fires.yr = fires[fires$YEAR == yrs[i],]#Select YROI
  fires.old = fires[fires$YEAR<yrs[i]&fires$YEAR>=(yrs[i]-window),]#Select prior 20 yrs
  print(Sys.time())
  fires.yr = rasterize(fires.yr,template)#Rasterize with 1ha pixel
  fires.old = rasterize(fires.old,template)
  comb = fires.yr+fires.old #Identify overlaps
  writeRaster(comb,paste(#Your path here
    "/Data/Yearly_Reburned_Area/Tiffs/Reburns_",window,"_",yrs[i],".tif",sep=""),overwrite=T)
  print(paste("Completed Year ",yrs[i],sep=""))
  i=i+1
}
