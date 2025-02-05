#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Feb 5, 2024
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

library(stars)
setwd(#Your path here
  "/Data/Yearly_Reburned_Area/Tiffs")
rasts = list.files(pattern = "\\.tif$") #List yearly reburn rasters
yrs = seq(1982,2023,by=1)#List of years for analysis
for(i in 1:length(rasts)){
  rb.ras = read_stars(rasts[i])#Read in
  rb.p2 = st_as_sf(rb.ras,as_points = FALSE,merge=T)#Convert to Poly using st
  rb.p2$YEAR = yrs[i]#Append year
  rm(rb.ras)
  gc()
  st_write(rb.p2, paste(#Your path here
    "/Data/Yearly_Reburned_Area/Polys/Reburns_10_",yrs[i],".shp",sep=""), delete_layer = TRUE)
  print(paste("stars Polygonized Year ",yrs[i],sep=""))
  i=i+1
}
rm(list=ls())
gc()

library(terra)
setwd(#Your path here
  "/Data/Yearly_Reburned_Area/Polys")
polys = list.files(pattern = "\\.shp$") #List yearly reburn polys
#This is inefficent, sorry, but we are now going to change this to terra polygons so that it is easier to calculate area, 
#and to make them separate patches for each year, and combine into a single polygon of the time series
yrs = seq(1982,2023,by=1)
for(j in 1:length(polys)){
  rb.pol = vect(polys[j])#Read in
  rb.pol = disagg(rb.pol)
  rb.pol$AREA_Ha = expanse(rb.pol,unit="ha")
  if(j==1){
    reburns = rb.pol
  }else{
    reburns = rbind(reburns,rb.pol)
  }
  print(paste("terra Polygonized Year ",yrs[j],sep=""))
  j=j+1
}
reburns = reburns[,c(2,3)]
head(reburns)

writeVector(reburns, 
            #Your path here
            "/Data/Reburns_TimeSeries_Complete_19802023_10yr.shp", overwrite = TRUE)
rm(list=ls())
gc()
