#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Feb 5, 2024
library(terra)
library(ggplot2)
library(Hmisc)
library(dplyr)
options(scipen = 999)

#The CMI & GDD Grids can be downloaded from the NRCan ftp.
#https://ftp.maps.canada.ca/pub/nrcan_rncan/Climate-archives_Archives-climatologiques/

#Import CMI grid and Reburn vector, remove small reburns, reproject all to EA Projection
cmi.20 = rast(#Your path here
  "mly60arcsecond_1991-2020/cmi60_sum.tif")
t.crs = "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Lambert Azimuthal Equal Area\",\n        METHOD[\"Lambert Azimuthal Equal Area\",\n            ID[\"EPSG\",9820]],\n        PARAMETER[\"Latitude of natural origin\",45,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]]]"
rebs = vect(#Your path here
  "/Data/Reburns_TimeSeries_Complete_19802023_20yr.shp")#This requires the polygon reburn maps created with scripts 1 & 2
rebs = rebs[rebs$AREA_Ha>=4,]#Remove small reburns
cmi.20 = project(cmi.20,crs(t.crs),method="bilinear")#Resampling CMI
rebs = project(rebs,crs(t.crs));rm(t.crs)

rebs.old = rebs[rebs$YEAR<2023,]#Reburns prior to 2023
rebs.old = rebs.old[rebs.old$YEAR>=1992,]#Truncate to reliable time period (accurate maps begin 1972, 20 year overlaps)
rebs.old.df = as.data.frame(rebs.old)
rebs = rebs[rebs$YEAR==2023,]#2023 only

rebs.old.df$ID = paste(row.names(rebs.old.df),rebs.old.df$YEAR,sep="_")#Create unique IDs for each old reburn
rebs.old$ID = rebs.old.df$ID

rebs.old.noseed = buffer(rebs.old,-200) #Interior buffer of >200m from edge/seed source
writeVector(rebs.old.noseed,
            #Your path here
            "/Data/Reburns_TimeSeries_Complete_19922023_20yr_InnerBuf.shp",overwrite=T)
#Read in existing seed-limited historical reburn vector file, if you are re-running this
#I think I maybe had some errors previously when I worked with the file as-generated, instead of re-importing the finalized file
#rebs.old.noseed = vect("/Users/ewhitman/Documents/Projects/Reburn_Occurence_2023/Data/Reburns_TimeSeries_Complete_19922023_20yr_InnerBuf.shp")

rebs.old.noseed$AREA_Noseed_Ha = expanse(rebs.old.noseed,"ha")#Calculate core area of each reburn in hectares
rebs.old.noseed.df = as.data.frame(rebs.old.noseed)
rebs.old = merge(rebs.old,rebs.old.noseed.df)#Some reburns have no core area, but we still want to know the original size. Merge together.
rebs.old$AREA_Seed_Ha = rebs.old$AREA_Ha-rebs.old$AREA_Noseed_Ha#Calculate non-core area
rebs.old.df = as.data.frame(rebs.old)#Convert from vector to data frame

#Aggregate the historical reburn data into years, calculate averages and standard devs
oldreburns.core.mean = rebs.old.df %>%
  group_by(YEAR) %>%
  summarize(Core_Area_Ha = sum(AREA_Noseed_Ha))
mean(oldreburns.core.mean$Core_Area_Ha)#34411.81
sd(oldreburns.core.mean$Core_Area_Ha)#31908.92


rebs.23.df = as.data.frame(rebs)
rebs.23.df$ID = paste(row.names(rebs.23.df),rebs.23.df$YEAR,sep="_")#Create unique IDs for 2023 reburns
rebs$ID = rebs.23.df$ID

rebs.23.noseed = buffer(rebs,-200)#Interior buffer of >200m from edge/seed source for 2023
writeVector(rebs.23.noseed,
            #Your path here
            "/Data/Reburns2023_Corearea.shp",overwrite=T)
#Read in existing seed-limited historical reburn vector if desired, as above
#rebs.23.noseed = vect("/Users/ewhitman/Documents/Projects/Reburn_Occurence_2023/Data/Reburns2023_Corearea.shp")

rebs.23.noseed$AREA_Noseed_Ha = expanse(rebs.23.noseed,"ha")#Calculate core area of each reburn in hectares
rebs.23.noseed.df = as.data.frame(rebs.23.noseed)
rebs = merge(rebs,rebs.23.noseed.df)#Some reburns have no core area, but we still want to know the original size. Merge together.
rebs$AREA_Seed_Ha = rebs$AREA_Ha-rebs$AREA_Noseed_Ha#Calculate non-core area
rebs.23.df = as.data.frame(rebs)#Convert from vector to data frame

#Calculate total reburn core area in 2023
sum(rebs.23.df$AREA_Noseed_Ha)#301536.1
(301536.1-34411.81)/31908.92 #8x outlier, 8.7X larger than historical mean
301536.1/34411.81

#We're only interested in 2023 for the plots, Create data frame for representing what happened in 2023
rebs.df = as.data.frame(rebs)
rebs.df$Area_Log = log10(rebs.df$AREA_Ha)
#range(log10(rebs.df$AREA_Ha))

#Set breakpoints for log area bins for barplot
temp2 = c(1.00000, 1.39794, 2.00000, 2.39794, 3.00000, 3.39794, 4.00000, 4.39794, 5.00000)
coreareadat <- transform(rebs.df, bin = cut2(Area_Log, temp2))

library(plyr)
coreareasum = ddply(coreareadat, "bin", summarize, 
                  AREA_Ha_Noseed = sum(AREA_Noseed_Ha),
                  AREA_Ha_Seed = sum(AREA_Seed_Ha),
                  Total_Patches = length(bin))
levels(coreareasum$bin)

#head(coreareasum)#Lookin' good
coreareasum$AREA_Ha = coreareasum$AREA_Ha_Noseed+coreareasum$AREA_Ha_Seed
coreareasum$AREA_Ha = as.numeric(coreareasum$AREA_Ha)
sum(coreareasum$AREA_Ha)#Does this still add up? Needed in order to calculate the cumulative sum 

#Calculate cumulative sum of reburned area in each column, and add to data frame
CumSum2 = mutate(coreareasum,Value = cumsum(AREA_Ha)) 
#head(coreareasum)
coreareasum$CumSum =CumSum2$Value
coreareasum$CumSumPct = (coreareasum$CumSum/1052583)*100
coreareasum$NLabel = paste("n=",coreareasum$Total_Patches,sep="")#Create label for number of reburns in each size class


core.temp = coreareasum[,c(1:3)]
names(core.temp) = c("bin","No_Seed","Seed")


library(tidyr)
corearea.long = core.temp %>%
  pivot_longer(!bin, names_to = "Seed_Risk", values_to = "AREA_Ha")#Reformat for bar plot

#Plot barplot of 2023 total area reburned, core area of reburn, by reburn size classes
corearea.size = ggplot(corearea.long, aes(fill=Seed_Risk, y=AREA_Ha, x=bin)) + 
  geom_bar(position="stack", stat="identity", width=1,colour = "black")+
  geom_text(data=coreareasum,aes(y = AREA_Ha+2000,x=bin,label=NLabel,fill = NA), vjust=0) +
  geom_line(data=coreareasum,aes(y = CumSumPct*2490,x=bin,group=1,fill=NA),colour="#355C7D",size=1)+
  scale_fill_manual(values=c("#E3B6BF","#C2ADC5"),guide="none")+
  xlab("Patch size (ha)")+
  ylab("Reburned area (ha)")+  
  scale_y_continuous(expand = expansion(mult=c(0,0.05)),
                     limits = c(0, 250000),
                     breaks = c(0,50000,100000,150000,200000,250000),
                     labels = c("0","50,000","100,000","150,000",
                                "200,000","250,000"),
    # Features of the first axis
    name = "Reburned area (ha)",
    # Add a second axis and specify its features
    sec.axis = sec_axis( transform=~./2490, name="Cumulative area (%)")) +
  scale_x_discrete(labels=c(">1",">5",">10",">25",">100",">250",">1,000",">2,500",
                            ">10,000",">25,000",">100,000"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill='transparent'),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.x= element_blank(),
        axis.text.y= element_text(colour="black",
                                  size=11),
        axis.ticks.x= element_blank(),)
corearea.size#Pretty!

#Make a template for rasterizing climate data with a 1ha pixel size
template = rast(ext(rebs),nlyrs=1,resolution = 100)
crs(template) = crs(rebs)
rebs.rast = rasterize(rebs,template)#Rasterize 2023 reburns, with 1ha resolution
#reb.temp.df = as.data.frame(rebs.rast)
#reb.temp.df$AREA_Ha = 1
#reb.temp.df = reb.temp.df[complete.cases(reb.temp.df),]
#sum(reb.temp.df$AREA_Ha)#This method loses us 43 ha of reburned area somewhere, because of going to raster at this lower res. This is equivalent to 0.004% of area reburned in Canada in 2023

#Assign CMI (climate moisture index) values from 1991-2020 normals to each area reburned in 2023
cmi.reb.20 = resample(cmi.20,rebs.rast,method="bilinear")
cmi.reb.20 = terra::mask(cmi.reb.20,rebs.rast)#Reburns only
#plot(cmi.reb.20)

#To data frame, for calculating total area reburned within CMI range/bins
reb.cmi.df = as.data.frame(cmi.reb.20)
reb.cmi.df$AREA_Ha = 1#Equal Area projection, each pixel is 1ha
reb.cmi.df = reb.cmi.df[complete.cases(reb.cmi.df),]
names(reb.cmi.df)[1]="CMI_2020"
#range(reb.cmi.df$CMI_2020)

#Create bins for barplot of reburned area CMI
temp = seq(-16,89,by=8)
temp = temp[-(10:14)]

#Create data frame for CMI of 2023 reburn barplots
cmidat <- transform(reb.cmi.df, bin = cut2(CMI_2020, temp))
cmisums = ddply(cmidat, "bin", summarize, AREA_Ha = sum(AREA_Ha))
#sum(cmisums$AREA_Ha)
CumSum = mutate(cmisums,Value = cumsum(AREA_Ha)) 
#head(cmisums)#Cool beans.
cmisums$CumSum =CumSum$Value
cmisums$CumSumPct = (cmisums$CumSum/1052540)*100

cmisums$CMI_Risk = NA
cmisums[1:2,"CMI_Risk"]= "Dry"
cmisums[3:nrow(cmisums),"CMI_Risk"]= "Moist"

#Plot area reburned by 2020s CMI
area.cmi = ggplot(cmisums)+
  geom_bar(aes(y = AREA_Ha,x=bin,fill = CMI_Risk),stat="identity",width=1,colour="black")+
  geom_line(aes(y = CumSumPct*4900,x=bin,group=1),colour="#355C7D", size=1)+
  #geom_vline(xintercept = 2.5)+
  scale_fill_manual(values=c("#FFC3A1","#FFAAA1"),guide="none")+ 
  scale_y_continuous(expand = expansion(mult=c(0,0.05)),
                     limits = c(0,500000),
                     breaks = c(0,100000,200000,300000,400000,500000),
                     labels = c("0","100,000","200,000","300,000",
                                "400,000","500,000"),
    # Features of the first axis
  name = "Reburned area (ha)",
    # Add a second axis and specify its features
   sec.axis = sec_axis( transform=~./4900, name="Cumulative area (%)")) +
  scale_x_discrete(labels=c(">-16",">-8",">0",">8",">16",
                            ">24",">32",">40",">48"))+
  xlab("Climate Moisture Index (1991 - 2020)")+
  ylab("Reburned area (ha)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill='transparent'),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y= element_text(colour="black",
                                  size=11))
area.cmi#Cute!

#Read in Growing Degree Days raster
#See above for download link
t.20 = rast(#Your path here
  "/mly60arcsecond_1991-2020/sg60_10.tif")
t.20 = project(t.20,crs(cmi.20),method="bilinear")
#Assign GDD (Growing Degree Day) values from 1991-2020 normals to each area reburned in 2023
t.reb.20 = resample(t.20,rebs.rast,method="bilinear")
t.reb.20 = terra::mask(t.reb.20,rebs.rast)#Only reburns

#Create data frame of GDD5 of in reburns
reb.gdd.df = as.data.frame(t.reb.20)
reb.gdd.df$AREA_Ha = 1 #Each pixel is 1 ha
reb.gdd.df = reb.gdd.df[complete.cases(reb.gdd.df),]
names(reb.gdd.df)[1]="DD5_2020"
#range(reb.gdd.df$DD5_2020)

#Breaks for barplot
temp = seq(200, 1600,by=200)

#Create data frame for GDD of area reburned analysis
gdd.bindat <- transform(reb.gdd.df, bin = cut2(DD5_2020, temp))
gdd.binsums = ddply(gdd.bindat, "bin", summarize, AREA_Ha = sum(AREA_Ha))
head(gdd.binsums)
#simplify by removing bins with extremely few pixels
gdd.binsums[2,2]=2829+11
gdd.binsums=gdd.binsums[-1,]
gdd.binsums[6,2]=1950+4
gdd.binsums=gdd.binsums[-7,]
sum(gdd.binsums$AREA_Ha)
gdd.binsums=gdd.binsums[-6,]

CumSum = mutate(gdd.binsums,Value = cumsum(AREA_Ha)) 
head(gdd.binsums)
gdd.binsums$CumSum =CumSum$Value
gdd.binsums$CumSumPct = (gdd.binsums$CumSum/1052540)*100

gdd.binsums$DD_Risk = NA
gdd.binsums[1:2,"DD_Risk"]= "Cold"
gdd.binsums[3:nrow(gdd.binsums),"DD_Risk"]= "Warm"

#Create plot of area reburned in 2023 by GDD of reburn
area.gdd = ggplot(gdd.binsums)+
  geom_bar(aes(y = AREA_Ha,x=bin,fill = DD_Risk),stat="identity",width=1,colour="black")+
  geom_line(aes(y = CumSumPct*6550,x=bin,group=1),colour="#355C7D",size=1)+
  #geom_vline(xintercept = 3.5)+
  scale_fill_manual(values=c("#B2D2EB","#7BADD1"),guide="none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.05)),
                     limits = c(0,660000),
                     breaks = c(0, 200000,400000,600000),
                     labels = c("0","200,000","400,000","600,000"),
    # Features of the first axis
    name = "Reburned area (ha)",
    # Add a second axis and specify its features
    sec.axis = sec_axis( transform=~./6550, name="Cumulative area (%)")) +
  scale_x_discrete(labels=c("<600","<800","<1,000","<1,200","<1,400",
                            "<1,600"))+
  xlab("Growing degree days (1991 - 2020)")+
  ylab("Reburned area (ha)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.background = element_rect(fill='transparent'),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.x= element_blank(),
        axis.ticks.x= element_blank(),
        axis.text.y= element_text(colour="black",
                                  size=11))
area.gdd#OoOoohhh!

#Create Figure of area reburned with limitations to regeneration
library(cowplot)
bottom_row <- cowplot::plot_grid(area.gdd, area.cmi, labels = c('b', 'c'), label_size = 14)

clim.rebs = cowplot::plot_grid(corearea.size, bottom_row,
                               labels = c("a",""),label_size = 14,
                               ncol = 1)
clim.rebs

#Save figure
png(#Your path here
  filename = "Fig4_Size_CMI_DD5_ofReburns_20yr.png",res=300,
    units = "cm", width = 25, height = 15.2)
clim.rebs
dev.off()

sum(rebs.23.noseed.df$AREA_Noseed_Ha)#301536.1 ha of reburn seed-limited
cmi.reb.20[cmi.reb.20>=0]=0
cmi.reb.20[cmi.reb.20<0]=1
global(cmi.reb.20,fun="sum",na.rm=T)#497,297 ha of reburn moisture-limited
t.reb.20[t.reb.20<800]=1
t.reb.20[t.reb.20>=800]=0
global(t.reb.20,fun="sum",na.rm=T)#47,003 ha of reburn temperature-limited

#Stupid invalid geometries need fixing.
library(sf)
new=st_read(#Your path here
  "Reburns2023_Corearea.shp")
new = st_make_valid(new)
new.vect = vect(new)

#Convert core area to raster, to combine with cold reburns and dry reburns
rebs.23.noseed.ras = terra::rasterize(new.vect,template)
rebs.rast[is.na(rebs.23.noseed.ras)]=0

#Merge all three resilience risk types, as raster
reb.any.risk = rebs.rast+cmi.reb.20+t.reb.20
plot(reb.any.risk)
reb.any.risk.1 = reb.any.risk
reb.any.risk.1[reb.any.risk.1>1]=1 #Reclassify so that it is 1 (any risk) or 0

global(reb.any.risk.1,fun="sum",na.rm=T)#658433 ha of reburn with any limitation


#Save the output, if you want, but uncommenting this section
# reb.any.risk.1[reb.any.risk.1==0]=NA
# writeRaster(reb.any.risk,"Reburned_Area_vulnerable.tif",
#             overwrite=T)
# writeRaster(reb.any.risk.1,"Reburned_Area_vulnerable_All1.tif",
#             overwrite=T)
# reb.any.risk.1 = rast("Reburned_Area_vulnerable_All1.tif")
# reb.any.risk.1 = as.polygons(reb.any.risk.1)
# writeVector(reb.any.risk.1,"Reburned_Area_vulnerable_All1.shp",
#             overwrite=T)


