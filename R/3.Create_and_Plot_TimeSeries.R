#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Displaying a time series of short-interval wildfire overlaps, and comparing 2023 to past years
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Nov 25, 2025
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file
library(terra)
library(dplyr)
library(ggplot2)
library(cowplot)
library(boot)

#Import fire maps and reburn maps
fires = vect(#your path here
  'NBAC_1972_2024_20250506.shp')#The download link for the NBAC fire perimeter maps is found in script 1
reburns = vect(#Your path here
  "/Data/Reburns_TimeSeries_Complete_19802024_20yr.shp")#This requires the reburn polygon map created with scripts 1 & 2
reburns$ID = paste(1:nrow(reburns),"_",reburns$YEAR,sep="")

#Equal area projection with m as horizontal unit
t.crs = "PROJCRS[\"unknown\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"Lambert Azimuthal Equal Area\",\n        METHOD[\"Lambert Azimuthal Equal Area\",\n            ID[\"EPSG\",9820]],\n        PARAMETER[\"Latitude of natural origin\",45,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-100,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]]]"
reburns = project(reburns,t.crs)
reburns$AREA_Ha = expanse(reburns,"ha")#Recalculating area in Equal Area Projection for consistency with tablular data later
reburns = reburns[reburns$AREA_Ha>4,] #Get rid of small patches (e.g. a single pixel)
sum(reburns[reburns$YEAR==2023,"AREA_Ha"])


#Categorical/factor year of reburning, separating out 2023 from past
reburns$CAT_YEAR = reburns$YEAR
reburns$CAT_YEAR = as.character(reburns$CAT_YEAR)
reburns[reburns$CAT_YEAR!="2023","CAT_YEAR"]="Reference"
reburns$CAT_YEAR = as.factor(reburns$CAT_YEAR)
#Truncate to reliable time period (accurate maps begin 1972, 20 year overlaps)
reburns = reburns[reburns$YEAR>=1992,]

#convert to data frame for time series plots
reburns.df = as.data.frame(reburns)
reburns.old = reburns.df[reburns.df$YEAR<2023,]
reburns.2023 = reburns.df[reburns.df$YEAR==2023,]

library(SiMRiv)
#Create bins for distribution of reburn patch sizes when logged
temp =  binCounts(log10(reburns.old$AREA_Ha), c(0.60,5), 45, log = FALSE)
temp = t(temp)
temp = seq(0.60,5,by=0.1)

#Calculate average number of patches by size from historical reburn data
yrs = seq(1992,2022,by=1)
for(i in 1:length(yrs)){
  yrdat = reburns.old[reburns.old$YEAR==yrs[i],]
  binsout = binCounts(log10(yrdat$AREA_Ha), c(0.60,5), 45, log = FALSE)
  binsout = as.data.frame(binsout)
  names(binsout)[1]="Count"
  binsout$Bin = temp
  binsout$YEAR = yrs[i]
  if(i==1){
    yearly.bins = binsout
  }else{
    yearly.bins=rbind(yearly.bins,binsout)
  }
  i=i+1
}

#head(yearly.bins)
#Calculate mean number of reburns by patch size class, from historical data
oldreburns.bins.mean = yearly.bins %>%
  group_by(Bin) %>%
  summarize(Count = mean(Count))

#Repeat for 2023, but no means (counts)
reburns2023.bins = binCounts(log10(reburns.2023$AREA_Ha), c(0.60,5), 45, log = FALSE)
reburns2023.bins = as.data.frame(reburns2023.bins)

#Prepare data frame for plot by creating factor and combining 2023 and historical data
names(reburns2023.bins)[1]="Count"
reburns2023.bins$Bin = temp
reburns2023.bins$YEAR = "2023"
oldreburns.bins.mean$YEAR = "Mean (1992 - 2022)"
reburns.binned = rbind(reburns2023.bins,oldreburns.bins.mean)
reburns.binned$YEAR = as.factor(reburns.binned$YEAR)

mean(reburns.old$AREA_Ha)#129.5 #Mean reburn patch size all historical data
#log10(129.5) #Convert to log
#Bootstrapped mean is the same

# Custom function mean
meanfun <- function(data, i){
  d <- data[i, ]
  return(mean(d))   
}
set.seed(42)
mean.patchsize <- boot(reburns.old[, "AREA_Ha", drop = FALSE], statistic=meanfun, R=1000)
mean.patchsize
boot.ci(mean.patchsize)


mean(reburns.2023$AREA_Ha)#161.6 #Mean reburn patch size in 2023
#log10(161.6)


max(reburns.2023$AREA_Ha)#120555 #Largest single patch
sum(reburns.2023$AREA_Ha)#1053249 #Total area reburned
120555/1053249

#Plot distributions of 2023 and mean historical reburned areas patch size
patchsize.dist<-ggplot(reburns.binned, aes(y=Count, x=Bin,group=YEAR,fill=YEAR,colour=YEAR)) +
  geom_density(alpha=0.4,stat = "identity")+
  scale_color_manual(values=c("#FFAAA1","#7BADD1"))+
  scale_fill_manual(values=c("#FFAAA1","#7BADD1"))+
  labs(y = "Number of reburn patches",x = "Size of reburned patch (ha)",fill="Time period",colour="Time period")+
  scale_x_continuous(expand = expansion(mult=c(0,0.05)),breaks=c(0,1,2,3,4,5),labels = c("0","10","100","1,000","10,000","100,000"))+
  geom_vline(xintercept = 2.11227,colour = "#7BADD1", linewidth = 0.7)+
  geom_vline(xintercept = 2.208441, colour = "#FFAAA1", linewidth = 0.7)+
  scale_y_continuous(expand = c(0,0), limits = c(0, 1000),labels=c("0","250","500","750","1,000"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(colour = "black"),
        legend.position = c(0.7, 0.88),
        legend.background = element_rect(fill='transparent',),
        legend.title=element_text(color="black",size=14),
        legend.text=element_text(color="black",size=12),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=13, 
                                    color="black"),
        axis.title.y = element_text(size=13, 
                                    color="black"),
        axis.text.x= element_text(colour="black",
                                  size=11),
        axis.text.y= element_text(colour="black",
                                  size=11))
patchsize.dist#You're not the boss of ME ggplot. I'll use deprecated terms if I want to!

#Calculate annual reburn summary statistics for each year 1992 - 2023
reburns.ann = reburns.df %>%
  group_by(YEAR) %>%
  summarize(Area_Reb_Ha = sum(AREA_Ha),
            N_Reb_Patches = n(),
            Biggest_Reb_Patch_Ha = max(AREA_Ha),
            Mean_Reb_Patch_Ha = mean(AREA_Ha),
            Median_Reb_Patch_Ha = median(AREA_Ha),
            P95_Reb_Patch_Ha = quantile(AREA_Ha,0.95))
#reburns.ann = reburns.ann[-33,]#getting rid of 2024 if included

#What about quite big patches (>100 ha in size)?
largepatches = reburns.df[reburns.df$AREA_Ha>100,]
reburns.lpcount = largepatches %>%
  group_by(YEAR) %>%
  summarize(N_Patches_100ha = n())
options(scipen=999)
#reburns.lpcount = reburns.lpcount[-33,]#getting rid of 2024 if included

#Set midpoint for gradient as 2008
mid = 2008 #2008 was mid.
reburns.ann = as.data.frame(reburns.ann)
mean(reburns.ann[1:31,"N_Reb_Patches"],na.rm=T)#Mean numpatches historically 797.5161
mean(reburns.ann[1:31,"Area_Reb_Ha"],na.rm=T)#Mean area reburned historically 103263.2

wilcox.test(reburns.old$AREA_Ha,reburns.2023$AREA_Ha)#Compare mean patch sizes. 2023 significantly > than historical

#Plot relationship between number of reburn patches and total area reburned for each year in time series
patchsize.num<-ggplot(reburns.ann,aes(x=N_Reb_Patches,y=Area_Reb_Ha))+
  geom_point(aes(colour=YEAR))+
  scale_color_gradient2(midpoint = mid, low = "#7BADD1", mid = "#C2ADC5",
                          high = "#FFAAA1")+
  annotate(geom = "point",x=797.5161,y=103263.2,colour="#355C7D",size = 3,shape = 17)+#Annotate historical mean values
  labs(x = "Number of reburn patches",y = "Total area reburned (ha)",colour=NULL)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1200000),
                     labels=c("0","300,000","600,000","900,000","1,200,000"))+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 7500),
                     labels=c("0","2,000","4,000","6,000",""))+
  annotate(geom="text", x=5800, y=1080000, 
             label="2023", color="black")+
  annotate(geom="text", x=2100, y=440000, 
           label="1995", color="black")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(colour = "black"),
        legend.position = c(0.19, 0.78),
        legend.background = element_rect(fill='transparent'),
        legend.title=element_text(color="black",size=14),
        legend.text=element_text(color="black",size=12),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=13, 
                                    color="black"),
        axis.title.y = element_text(size=13, 
                                    color="black"),
        axis.text.x= element_text(colour="black",
                                  size=11),
        axis.text.y= element_text(colour="black",
                                  size=11))
patchsize.num#Oooh la la

#Calculate area burned by fire (not exclusively reburns) by year
fires.df = as.data.frame(fires)
allfires.ann = fires.df %>%
  group_by(YEAR) %>%
  summarize(N_Fires = n(),
            Area_Burned_Total_Ha = sum(POLY_HA))
allfires.ann = as.data.frame(allfires.ann)


#Create time series of reburn data
dat.ts = as.data.frame(reburns.ann)
mean(dat.ts[1:31,"Area_Reb_Ha"])#103263.2 Average area reburned prior to 2023
max(dat.ts[1:31,"Area_Reb_Ha"])#410097.5 previous max area reburned (in 1995)
max(dat.ts[1:32,"Area_Reb_Ha"])#1053249 Area reburned 2023
1053249/103263.2#10.2 x larger
1053249/410097.5#2.6 x larger

set.seed(42)
mean.arb <- boot(dat.ts[1:31, "Area_Reb_Ha", drop = FALSE], statistic=meanfun, R=1000)#Excluding 2023
mean.arb
boot.ci(mean.arb)

#Not used in analysis/manuscript but can calculate how extreme of an outlier 2023 was
# sd(dat.ts[1:31,"Area_Reb_Ha"]) #historical SD 88434.83
# (1053249-103263.2)/88434.83 #area reburned in 2023 is a 10.7 X outlier


mean(dat.ts[1:31,"N_Reb_Patches"])#797.5161 average number of reburn patches prior to 2023
max(dat.ts[1:31,"N_Reb_Patches"])#2751 previous max number, also 1995
max(dat.ts[1:32,"N_Reb_Patches"])#6519 number of reburn patches 2023

# sd(dat.ts[1:31,"N_Reb_Patches"]) #sd 564.6706
# (6519-797.5161)/564.6706 #number of reburn patches in 2023 10.1 X outlier

mean(dat.ts[1:31,"Mean_Reb_Patch_Ha"])#125.3191 historical average reburn patch size#Mean of mean, vs. mean of all
max(dat.ts[1:31,"Mean_Reb_Patch_Ha"])#250.9822 historical previous max number #2001! Twist!
mean(dat.ts[32,"Mean_Reb_Patch_Ha"]) #161.5661 2023 avg Did not need put this is in mean() but oh well

# sd(dat.ts[1:31,"Mean_Reb_Patch_Ha"]) #sd 46.70467
# (161.5661-125.31916)/46.70467 #0.78 X outlier. Not that impressive.

mean(dat.ts[1:31,"Biggest_Reb_Patch_Ha"])#14443.96 average biggest reburn patch size
max(dat.ts[1:31,"Biggest_Reb_Patch_Ha"])#50900.77 previous max number #1993
max(dat.ts[1:32,"Biggest_Reb_Patch_Ha"])#120555 max patch size 2023

# sd(dat.ts[1:31,"Biggest_Reb_Patch_Ha"]) #sd  13543.86
# (120555-14443.96)/ 13543.86 #2023 largest patch is 7.8 X outlier

#More stats about reburns. Creating data frame with total AB and largest single patch for time series plot
dat.ts = merge(allfires.ann,reburns.ann)
dat.ts = merge(dat.ts,reburns.lpcount)
dat.ts$Pct_Biggest = (dat.ts$Biggest_Reb_Patch_Ha/dat.ts$Area_Reb_Ha)*100
dat.ts$Pct_Reburn = (dat.ts$Area_Reb_Ha/dat.ts$Area_Burned_Total_Ha)*100
dat.ts$Area_Reburns_MinusLargest = dat.ts$Area_Reb_Ha-dat.ts$Biggest_Reb_Patch_Ha
dat.ts$Prop_Largepatch = dat.ts$N_Patches_100ha / dat.ts$N_Reb_Patches

dat.reburn.long = dat.ts[,c("YEAR","Area_Reburns_MinusLargest")]
dat.reburn.long$ReburnType = "All_Reburn"
names(dat.reburn.long)[2]="Area_Reburned_Ha"
temp = dat.ts[,c("YEAR","Biggest_Reb_Patch_Ha")]
temp$ReburnType = "Largest_Patch"
names(temp)[2]="Area_Reburned_Ha"

dat.reburn.long = rbind(dat.reburn.long,temp);rm(temp)
options(scipen = 999)

mean(dat.ts[1:31,"Area_Reb_Ha"]) #historical mean total area reburned 103263.2
#Plot value as horitzonal line on barplot

#Plot time series of total area reburned and proportion contributed by single largest patch
ReburnArea.Ann = ggplot(dat.reburn.long, aes(y=Area_Reburned_Ha, x=YEAR,fill=ReburnType)) +
  geom_bar(position="stack", stat="identity",colour="black")+
  scale_fill_manual(values=c("#C2ADC5","#E3B6BF"))+
  geom_hline(yintercept = 103263.2,linewidth = 0.7,linetype="dashed",colour="#355C7D")+
  labs(y = "Area reburned (ha)",x = "Year")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA),labels=c("0","250,000","500,000","750,000","1,000,000"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=13, 
                                    color="black"),
        axis.title.y = element_text(size=13, 
                                    color="black"),
        axis.text.x= element_text(colour="black",
                                  size=11),
        axis.text.y= element_text(colour="black",
                                  size=11))
ReburnArea.Ann

bottom_row = cowplot::plot_grid(patchsize.dist, patchsize.num, labels = c('b', 'c'), label_size = 15)

reburns.alldat = cowplot::plot_grid(ReburnArea.Ann, bottom_row, labels = c('a', ''), 
                                    label_size = 15, ncol = 1)
reburns.alldat

#Save figure
png(#Your path here
  filename = "/Users/ewhitman/Documents/Projects/Reburn_Occurence_2023/Revision/Figures/Annual_Area_Reburn_BarChart_Sizedist_Npatch_20yr.png",res=300,
    units = "cm", width = 21, height = 16)
reburns.alldat
dev.off()

library(trend)
mk.test(dat.ts[1:32,"Area_Reb_Ha"],alternative = "greater") #p = 0.15

#Spearman's Correlation tests between annual fire and reburn data
cor.test(log10(dat.ts$Area_Reb_Ha),log10(dat.ts$Area_Burned_Total_Ha),
         alternative = c("two.sided"),
         method = c("spearman"))

ab.abreb = ggplot(dat.ts,aes(x = log10(Area_Burned_Total_Ha), y = log10(Area_Reb_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=7.164354, y=6.022531, color="#FFAAA1", size = 1.8)+
  ylab("log(Area reburned (ha))")+
  xlab("log(Total area burned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
ab.abreb

cor.test(log10(dat.ts$Area_Reb_Ha),dat.ts$Pct_Reburn,
         alternative = c("two.sided"),
         method = c("spearman"))

pctreb.abreb = ggplot(dat.ts,aes(y = Pct_Reburn, x = log10(Area_Reb_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=6.022531, y=7.214015 , color="#FFAAA1", size = 1.8)+
  xlab("log(Area reburned (ha))")+
  ylab("Area burned as reburn (%)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
pctreb.abreb

cor.test(log10(dat.ts$Area_Burned_Total_Ha),dat.ts$Pct_Reburn*100,
         alternative = c("two.sided"),
         method = c("spearman"))

pctreb.ab = ggplot(dat.ts,aes(y = Pct_Reburn, x = log10(Area_Burned_Total_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=7.164354, y=7.214015, color="#FFAAA1", size = 1.8)+
  ylab("Area burned as reburn (%)")+
  xlab("log(Total area burned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
pctreb.ab

cor.test(log10(dat.ts$Area_Burned_Total_Ha),log10(dat.ts$Biggest_Reb_Patch_Ha),
         alternative = c("two.sided"),
         method = c("spearman"))

maxpatch.ab = ggplot(dat.ts,aes(y = log10(Biggest_Reb_Patch_Ha), x = log10(Area_Burned_Total_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=7.164354, y=5.081185, color="#FFAAA1", size = 1.8)+
  ylab("log(Largest reburn patch (ha))")+
  xlab("log(Total area burned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
maxpatch.ab

cor.test(log10(dat.ts$Area_Reb_Ha),log10(dat.ts$Biggest_Reb_Patch_Ha),
         alternative = c("two.sided"),
         method = c("spearman"))

maxpatch.abreb = ggplot(dat.ts,aes(y = log10(Biggest_Reb_Patch_Ha), x = log10(Area_Reb_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=6.022531, y=5.081185, color="#FFAAA1", size = 1.8)+
  ylab("log(Largest reburn patch (ha))")+
  xlab("log(Area reburned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
maxpatch.abreb

cor.test(log10(dat.ts$Area_Burned_Total_Ha),dat.ts$N_Reb_Patches,
         alternative = c("two.sided"),
         method = c("spearman"))

npatch.ab = ggplot(dat.ts,aes(y = N_Reb_Patches, x = log10(Area_Burned_Total_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=7.164354, y=6519, color="#FFAAA1", size = 1.8)+
  ylab("Number of patches")+
  xlab("log(Total area burned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
npatch.ab

cor.test(log10(dat.ts$Area_Reb_Ha),dat.ts$N_Reb_Patches,
         alternative = c("two.sided"),
         method = c("spearman"))

npatch.abreb = ggplot(dat.ts,aes(y = N_Reb_Patches, x = log10(Area_Reb_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=6.022531, y=6519, color="#FFAAA1", size = 1.8)+
  ylab("Number of patches")+
  xlab("log(Area reburned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
npatch.abreb

cor.test(log10(dat.ts$Area_Reb_Ha),dat.ts$N_Patches_100ha,
         alternative = c("two.sided"),
         method = c("spearman"))

lgpatch.abreb = ggplot(dat.ts,aes(y = N_Patches_100ha, x = log10(Area_Reb_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=6.022531, y=820, color="#FFAAA1", size = 1.8)+
  ylab("Number of large patches")+
  xlab("log(Area reburned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
lgpatch.abreb

cor.test(log10(dat.ts$Area_Burned_Total_Ha),dat.ts$N_Patches_100ha,
         alternative = c("two.sided"),
         method = c("spearman"))

lgpatch.ab = ggplot(dat.ts,aes(y = N_Patches_100ha, x = log10(Area_Burned_Total_Ha)))+
  geom_point(colour="#355C7D")+
  annotate(geom="point", x=7.164354, y=820, color="#FFAAA1", size = 1.8)+
  ylab("Number of large patches")+
  xlab("log(Total area burned (ha))")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"))
lgpatch.ab

corplots = cowplot::plot_grid(ab.abreb, pctreb.ab, pctreb.abreb,
                              npatch.ab,npatch.abreb,lgpatch.ab,
                              lgpatch.abreb,maxpatch.ab, maxpatch.abreb, 
                              labels = c('a','b','c',
                                         'd','e','f',
                                         'g','h','i'),
                              nrow=3)
corplots

#Save figure
png(#Your path here
  filename = "/Users/ewhitman/Documents/Projects/Reburn_Occurence_2023/Revision/Figures/S2_Correlations_Reburn_Areaburned_Pctreburn.png",res=300,
    units = "cm", width = 25, height = 18)
corplots
dev.off()
