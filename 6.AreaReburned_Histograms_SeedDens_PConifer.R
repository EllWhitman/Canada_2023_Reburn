#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Analysis of climate and core area within 2023 short-interval reburned areas
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Nov 25, 2025
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

rm(list=ls())
library(terra)
library(ggplot2)
library(glmmTMB)
library(Hmisc)
library(plyr)
library(cutpointr)
options(scipen = 999)

#Read in reburns polygon
reburns = vect(#Your path here
  "/Data/Reburns_TimeSeries_Complete_19802024_20yr.shp")
reburns$ID = paste(1:nrow(reburns),"_",reburns$YEAR,sep="")#Make unique IDs for each patch
t.crs = "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Lambert Azimuthal Equal Area\",\n        METHOD[\"Lambert Azimuthal Equal Area\",\n            ID[\"EPSG\",9820]],\n        PARAMETER[\"Latitude of natural origin\",45,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]]]"
reburns = project(reburns,t.crs)#Make same projection as rest of data (equal area)
reburns$AREA_Ha = expanse(reburns,"ha")#Calculate area in correct projection
reburns = reburns[reburns$YEAR==2023,]#2023 only

corearea = buffer(reburns,-300)#Inner buffer of >300m from edge of reburn (core area)
#Core area thresholds identified in modelling script
corearea$Area_noseed_Ha = expanse(corearea,"ha")#area of corearea

#Convert to data frame and merge by unique ID so that core area is matched to the same reburn patch
corearea = as.data.frame(corearea)
reburns = as.data.frame(reburns)
reburns = merge(reburns,corearea)
reburns$Area_seed_Ha = reburns$AREA_Ha-reburns$Area_noseed_Ha
reburns.2023 = reburns[reburns$AREA_Ha>=4,]#Limit to larger reburn patches, as throughout

#Core area of reburns lacking seed sources
sum(reburns.2023$Area_noseed_Ha)#178244
sum(reburns.2023$AREA_Ha)#1053249
#Core area of large reburns
sum(reburns.2023[reburns.2023$AREA_Ha>=100,"Area_noseed_Ha"])#178209.5
sum(reburns.2023[reburns.2023$AREA_Ha>=14200,"Area_noseed_Ha"])#91123.94
nrow(reburns.2023[reburns.2023$AREA_Ha>=100,])
nrow(reburns.2023)
# 178244/1053249
# 178209.5/178244
# 820/6519
# 91123.94/178244

rebs.df = reburns.2023#Convert to data frame
#Calculate core area by reburn patch size and display as barplot
rebs.df$Area_Log = log10(rebs.df$AREA_Ha)

#Set breakpoints for log area bins for barplot, and sum area within bins
temp2 = c(1.00000, 1.39794, 2.00000, 2.39794, 3.00000, 3.39794, 4.00000, 4.39794, 5.00000)
coreareadat <- transform(rebs.df, bin = cut2(Area_Log, temp2))

library(plyr)
coreareasum = ddply(coreareadat, "bin", summarize, 
                    AREA_Ha_Noseed = sum(Area_noseed_Ha),
                    AREA_Ha_Seed = sum(Area_seed_Ha),
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
coreareasum$CumSumPct = (coreareasum$CumSum/1053249)*100#cumsum is by 1053249 becasue these are the original polygons
coreareasum$NLabel = paste("n=",coreareasum$Total_Patches,sep="")#Create label for number of reburns in each size class

core.temp = coreareasum[,c(1:3)]
names(core.temp) = c("bin","No_Seed","Seed")

library(tidyr)
#Pivot to longer/reformat for a barplot
corearea.long = core.temp %>%
  pivot_longer(!bin, names_to = "Seed_Risk", values_to = "AREA_Ha")

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

#Tidy up
rm(core.temp,corearea,coreareadat,coreareasum,old.reburns,reburns,reburns.2023,temp2)

#Read in predicted post-fire trajectory data as .csv
climdat = read.csv(#Your path here
  "/Data/Regen_Trajectories_Data.csv")
climdat$Area_Ha = 1#Raster has a resolution where 1 cell = 1ha, therefore 1 row = 1 ha

#range(climdat$CMI)#-14.3 to 88.9
#Breakpoints for CMI bins
temp = seq(-16,89,by=8)
temp = temp[-(10:14)]#Truncate to remove long tails

#Create data frame for CMI of 2023 reburn barplots
cmidat <- transform(climdat, bin = cut2(CMI, temp))
#Sum of area of each trajectory by CMI bin
cmi.con = ddply(cmidat[cmidat$Trajectory=="Conifer",], "bin", summarize, AREA_Ha = sum(Area_Ha))
cmi.broad = ddply(cmidat[cmidat$Trajectory=="Broadleaf",], "bin", summarize, AREA_Ha = sum(Area_Ha))
cmi.sparse = ddply(cmidat[cmidat$Trajectory=="Sparse",], "bin", summarize, AREA_Ha = sum(Area_Ha))
cmisums2 = ddply(cmidat, "bin", summarize, AREA_Ha = sum(Area_Ha))
cmi.con$Trajectory = "Conifer"
cmi.broad$Trajectory = "Broadleaf"
cmi.sparse$Trajectory = "Sparse"
#Combine into one data frame
cmisums = rbind(cmi.con,cmi.broad)
cmisums = rbind(cmisums,cmi.sparse)

#sum(cmisums$AREA_Ha)
CumSum = mutate(cmisums2,Value = cumsum(AREA_Ha)) 
#head(cmisums)#Cool beans.
cmisums2$CumSum =CumSum$Value
cmisums2$CumSumPct = (cmisums2$CumSum/1053242)*100 #cumsum is by 1053242 becasue we lose 7 ha converting to raster

#Plot area reburned by binned 2020s CMI
area.cmi = ggplot(cmisums)+
  geom_bar(aes(y = AREA_Ha,x=bin,fill = Trajectory),stat="identity",width=1,colour="black")+
  geom_line(data = cmisums2,aes(y = CumSumPct*4900,x=bin,group=1),colour="#355C7D", size=1)+
  scale_fill_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.05)),
                     limits = c(0,500000),
                     breaks = c(0,100000,200000,300000,400000,500000),
                     labels = c("0","100,000","200,000","300,000",
                                "400,000","500,000"),
                     # Features of the first axis
                     name = "",
                     # Add a second axis and specify its features
                     sec.axis = sec_axis( transform=~./4900, name="Cumulative area (%)")) +
  xlab("Climate Moisture Index (1991 - 2020)")+
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

#Not used in final analysis but if you want to calculate the proportion of each trajecotry by bin, you can.
# library(ggstats)
# cmi.p <- ggplot(cmisums) +
#   aes(x = bin, fill = Trajectory,weight=AREA_Ha,by=bin) +
#   geom_bar(position = "fill") +
#   geom_text(stat = "prop", position = position_fill(.5))+
#   scale_fill_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
#   xlab("Climate Moisture Index (1991 - 2020)")+
#   ylab("Proportion")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.background = element_rect(fill='transparent'),
#         text = element_text(colour = "black"),
#         axis.title.x = element_text(size=14, 
#                                     color="black"),
#         axis.title.y = element_text(size=14, 
#                                     color="black"),
#         axis.text.x= element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y= element_text(colour="black",
#                                   size=11))
# cmi.p


#Breaks for GDD barplot
temp = seq(200, 1600,by=100)

#Create data frame for analysis of GDD within short-interval reburned area
gdd.bindat <- transform(climdat, bin = cut2(GDD, temp))
#Sums of area by predicted vegetation trajectory
gdd.con = ddply(gdd.bindat[gdd.bindat$Trajectory=="Conifer",], "bin", summarize, AREA_Ha = sum(Area_Ha))
gdd.broad = ddply(gdd.bindat[gdd.bindat$Trajectory=="Broadleaf",], "bin", summarize, AREA_Ha = sum(Area_Ha))
gdd.sparse = ddply(gdd.bindat[gdd.bindat$Trajectory=="Sparse",], "bin", summarize, AREA_Ha = sum(Area_Ha))
gdd.binsums2 = ddply(gdd.bindat, "bin", summarize, AREA_Ha = sum(Area_Ha))
gdd.con$Trajectory = "Conifer"
gdd.broad$Trajectory = "Broadleaf"
gdd.sparse$Trajectory = "Sparse"

#simplify by removing bins with extremely few pixels
gdd.broad[2,"AREA_Ha"] = 4524
gdd.broad[9,"AREA_Ha"] = 1242
gdd.broad = gdd.broad[-c(1,10,11,12),]

gdd.con[4,"AREA_Ha"] = 3259
gdd.con[11,"AREA_Ha"] = 799
gdd.con = gdd.con[-c(1,2,3,12,13),]

gdd.sparse[3,"AREA_Ha"] = 3508
gdd.sparse[10,"bin"] = "[1300,1400)"
gdd.sparse=gdd.sparse[-c(1,2),]

gdd.binsums = rbind(gdd.con,gdd.broad)
gdd.binsums = rbind(gdd.binsums,gdd.sparse)
head(gdd.binsums)

gdd.binsums2[4,"AREA_Ha"] = 11291
gdd.binsums2[11,"AREA_Ha"] = 2059
gdd.binsums2=gdd.binsums2[-c(1:3,12:14),]

#Calculate cumulative sum by bin
CumSum = mutate(gdd.binsums2,Value = cumsum(AREA_Ha)) 
head(gdd.binsums2)
gdd.binsums2$CumSum =CumSum$Value
gdd.binsums2$CumSumPct = (gdd.binsums2$CumSum/1053242)*100


#Create plot of area reburned in 2023 by GDD of reburn
area.gdd = ggplot(gdd.binsums)+
  geom_bar(aes(y = AREA_Ha,x=bin,fill = Trajectory),stat="identity",width=1,colour="black")+
  geom_line(data = gdd.binsums2,aes(y = CumSumPct*4950,x=bin,group=1),colour="#355C7D",size=1)+
  scale_fill_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
  scale_y_continuous(expand = expansion(mult=c(0,0.05)),
                     limits = c(0,500000),
                     breaks = c(0, 100000,200000,300000,400000,500000),
                     labels = c("0","100,000","200,000","300,000","400,000","500,000"),
                     # Features of the first axis
                     name = "Reburned area (ha)",
                     # Add a second axis and specify its features
                     sec.axis = sec_axis( transform=~./4950, name="")) +
  xlab("Growing degree days (1991 - 2020)")+
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

# gdd.p <- ggplot(gdd.binsums) +
#   aes(x = bin, fill = Trajectory,weight=AREA_Ha,by=bin) +
#   geom_bar(position = "fill") +
#   geom_text(stat = "prop", position = position_fill(.5))+
#   scale_fill_manual(values=c('#C2ADC5','#7BADD1','#FFAAA1'),guide="none")+
#   xlab("Growing degree days (1991 - 2020)")+
#   ylab("Proportion")+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"),
#         legend.background = element_rect(fill='transparent'),
#         text = element_text(colour = "black"),
#         axis.title.x = element_text(size=14, 
#                                     color="black"),
#         axis.title.y = element_text(size=14, 
#                                     color="black"),
#         axis.text.x= element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.y= element_text(colour="black",
#                                   size=11))
# gdd.p

library(cowplot)
bottom_row <- cowplot::plot_grid(area.gdd, area.cmi, labels = c('b', 'c'), label_size = 14)

clim.rebs = cowplot::plot_grid(corearea.size, bottom_row,
                               labels = c("a",""),label_size = 14,
                               ncol = 1)
clim.rebs


png(filename = #Your path here
      "/Figures/Fig4_Size_CMI_DD5_ofReburns_20yr.png",res=300,
    units = "cm", width = 25, height = 17)
clim.rebs
dev.off()

#Summaries of total area (nrows, reminder 1 row = 1 ha) of reburn across different
nrow(climdat[climdat$CMI<=0,])#Reburned area that is dry
nrow(climdat[climdat$GDD<=800,])#Reburned area that is cool
nrow(climdat[climdat$Trajectory=="Sparse"&climdat$Distance_Reburnedge_m>300,])#Total area sparse that is >300m from the edge
range(climdat[climdat$Trajectory == "Conifer","CMI"])#Range of CMI with predicted conifer trajectory
range(climdat[climdat$Trajectory == "Conifer","GDD"])#Range of GDD with predicted conifer trajectory
range(climdat[climdat$Trajectory == "Conifer","Distance_Reburnedge_m"])#Range of distance from edge with predicted conifer trajectory
range(climdat$GDD)#Range of all GDD
table(climdat$Trajectory)#Broadleaf = 683696, 0.65; conifer = 162207, 0.15; sparse = 207339 0.20
table(climdat[climdat$CMI<=0|climdat$GDD<=800&climdat$Distance_Reburnedge_m<=300,"Trajectory"])# Broadleaf = 403173, 0.755; Conifer = 8029,0.015; Sparse = 122638 0.23 
table(climdat[climdat$CMI<=0|climdat$GDD<=800&climdat$Distance_Reburnedge_m>300,"Trajectory"])# Broadleaf = 377026, 0.74; Conifer = 105,0.0; Sparse = 129527 0.26
table(climdat[climdat$CMI<=0|climdat$GDD<=800,"Trajectory"])#Broadleaf = 393882     Conifer = 8018    Sparse = 141286

nrow(climdat[(climdat$CMI<=0|climdat$GDD<=800)&climdat$Total_Stems_ha<=100,])#Area very sparse, in cool and dry climates
nrow(climdat[climdat$Total_Stems_ha<=100,])#Total area very sparse
nrow(climdat[climdat$Total_Stems_ha<=100&climdat$Distance_Reburnedge_m<=300,])#Area very sparse far from edge

#Area of each trajectory in dry and cool climates
clim.obs = c(393882+8018,141286)
#Proportion of each trajectory in full dataset
clim.exp = c(0.65+0.15,0.2)

chisq.test(x = clim.obs,p = clim.exp)#P<0.001,  sparse trajectories are more likely below a CMI of 0 or GDD of 800