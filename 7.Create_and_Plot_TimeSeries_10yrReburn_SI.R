#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Nov 25, 2025
library(terra)
library(dplyr)
library(ggplot2)
library(cowplot)
library(trend)

#Re-run scripts 1 & 2, but with a 10-year window between fires, in order to create the reburn time series
#file required in this first line of code.

reburns = vect(#Your path and filename here
  "Reburns_TimeSeries_Complete_19822023_10yr.shp")
reburns$ID = paste(1:nrow(reburns),"_",reburns$YEAR,sep="")

reburns$CAT_YEAR = reburns$YEAR
reburns$CAT_YEAR = as.character(reburns$CAT_YEAR)
reburns[reburns$CAT_YEAR!="2023","CAT_YEAR"]="Reference"
reburns$CAT_YEAR = as.factor(reburns$CAT_YEAR)
reburns = reburns[reburns$YEAR>=1982,]#Truncate to valid period (10-year overlap window, comprehensive maps begin 1972)

reburns.df = as.data.frame(reburns)
reburns.df = reburns.df[reburns.df$AREA_Ha>4,]
reburns.old = reburns.df[reburns.df$YEAR<2023,]#Divide into historical
reburns.2023 = reburns.df[reburns.df$YEAR==2023,]#And 2023

reburns.ann = reburns.df %>%#Calculate annual summary statistics
  group_by(YEAR) %>%
  summarize(Area_Reb_Ha = sum(AREA_Ha),
            N_Reb_Patches = n(),
            Biggest_Reb_Patch_Ha = max(AREA_Ha),
            Mean_Reb_Patch_Ha = mean(AREA_Ha),
            Median_Reb_Patch_Ha = median(AREA_Ha),
            P95_Reb_Patch_Ha = quantile(AREA_Ha,0.95))


sens.slope(reburns.ann$Area_Reb_Ha,conf.level=0.95)
mk.test(reburns.ann$Area_Reb_Ha,alternative = "greater")#Significant increase in area burned over time. B = 737.0747439379 p < 0.001
m=sens.slope(reburns.ann$Area_Reb_Ha,conf.level=0.95)
m=m$estimates
b=median(reburns.ann$Area_Reb_Ha)-m*median(reburns.ann$YEAR)#Calculate slope and intercept for representing Sen's slope
reburns.ann$Reburn10yr_AB_lm=reburns.ann$YEAR*m+b#Add linear model values for each year

ReburnArea.Ann = ggplot(reburns.ann, aes(y=Area_Reb_Ha, x=YEAR)) +
  geom_bar(position="stack", stat="identity",fill="#C2ADC5")+
  geom_line(aes(YEAR,Reburn10yr_AB_lm),linewidth = 0.7,colour="#355C7D")+
  labs(y = "Area reburned (ha)",x = "Year")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA),labels=c("0","100,000","200,000","300,000","400,000","500,000"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(1,1,1.5,1), "lines"),
        text = element_text(colour = "black"),
        axis.title.x = element_text(size=14, 
                                    color="black"),
        axis.title.y = element_text(size=14, 
                                    color="black"),
        axis.text.x= element_text(colour="black",
                                  size=11),
        axis.text.y= element_text(colour="black",
                                  size=11))
ReburnArea.Ann

#Save figure
png(filename = 
      #Your path here
      "FigS2_Annual_Area_Reburn_BarChart_10yr_SI.png",res=300,
    units = "cm", width = 16, height = 10)
ReburnArea.Ann
dev.off()

library(trend)
mk.test(reburns.ann$Area_Reb_Ha,alternative = "greater") #p = 0.015
mk.test(reburns.ann$N_Reb_Patches,alternative = "greater") #p = 0.001
