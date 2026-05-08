#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Model validation and transferrablity metrics. Model assessment.
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit May 8, 2026
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

rm(list=ls())
library(ggplot2)
library(stars)
library(glmmTMB)
library(ROSE)
library(blockCV)
library(performance)
library(Metrics)
library(pscl)
library(terra)

#Model Assessment and validation
#Calculate cross-validated model performance metrics and characterize performance
#Add spatial cross-validation
options(scipen = 999)
#Import field data and create response and predictor variables
#Same data processing as in all other scripts
seed.df = read.csv(#Your path here
  "/Data/Seedling_Data_Processed_ClimDist_March2326.csv")
seed.df = seed.df[seed.df$Total_seed_ha<400000,]
seed.df[seed.df$Interval=="Long","Distance_Reburnedge_m"]=0
seed.df$Interval=as.factor(seed.df$Interval)
seed.df$Broad_seed_ha = seed.df$Total_seed_ha - seed.df$Con_seed_ha
seed.df$Pct_Con = seed.df$Con_seed_ha/seed.df$Total_seed_ha*100
seed.df[seed.df$Total_seed_ha==0,"Pct_Con"]=0
seed.df$ConRegen = NA
seed.df[seed.df$Pct_Con>=61,"ConRegen"]=1 #Conifers dominate post-fire
seed.df[is.na(seed.df$ConRegen),"ConRegen"]=0 #Conifers do not dominate post-fire
seed.df$Sparse = NA
seed.df[seed.df$Total_seed_ha<=1000,"Sparse"]=1#Sparse regen in the field
seed.df[is.na(seed.df$Sparse),"Sparse"]=0
seed.df=seed.df[,-19]


#How many tree seedlings can establish after fire?
pseed = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
                ziformula=~Distance_Reburnedge_m+GDD,
                family=poisson(),
                data=seed.df)
summary(pseed)
pseed2 = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
                 data=seed.df,
                 family=poisson())#Model with no ZI
pseed3 = glmmTMB(Total_seed_ha~GDD+CMI+Distance_Reburnedge_m+years_post,
                 data=seed.df,
                 family=poisson())#Model with no random effect
check_overdispersion(pseed)#All models are overdispersed except the ZIP model
check_overdispersion(pseed2)
check_overdispersion(pseed3)

AIC(logLik(pseed))
AIC(logLik(pseed2))
AIC(logLik(pseed3))#AIC of Pseed with random effect and ZI is lowest!

pR2(pseed)#McFadden R2 of both random effects models is comparable
pR2(pseed2)
pR2(pseed3)
rm(pseed2,pseed3)

#Model validation
#How well do the predicted stem densities do, relative to the actual?
mae(seed.df$Total_seed_ha,predict(pseed,type="response"))#25120.58 #Actual MAE
mae(seed.df$Total_seed_ha,rep(mean(seed.df$Total_seed_ha),171))#36042.38 #Naive baseline

rmsle(seed.df$Total_seed_ha,predict(pseed,type="response"))#2.56 #Actual RMSLE
rmsle(seed.df$Total_seed_ha,rep(mean(seed.df$Total_seed_ha),171))#2.9 #Naive baseline

mean(abs(log10((seed.df$Total_seed_ha+1)/(predict(pseed,type="response")+1))))#0.68 #Actual MALE
mean(abs(log10((seed.df$Total_seed_ha+1)/(rep(mean(seed.df$Total_seed_ha),171)+1))))#0.84 #Naive baseline

#The data is pretty imbalanced, so we will create a new dataset with random up and downsampling
#For modelling conifer dominance likelihood
set.seed(1234)
seed.over = ovun.sample(ConRegen~., data = seed.df, method = "both", N = 200)
table(seed.over$data$ConRegen)

#Model of likelihood of post-fire conifer dominance
pcon = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m,
               data=seed.over$data,
               family=binomial)#Simple model
pcon2 = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m+(1|project_id),
                data=seed.over$data,
                family=binomial)#Mixed effects model

AIC(logLik(pcon))
AIC(logLik(pcon2))#AIC is improved without the random effect. User the simpler model.
check_overdispersion(pcon)#Neither model is overdispersed
check_overdispersion(pcon2)

pR2(pcon)#Pr2s are identical
pR2(pcon2)
rm(pcon2)

#Predict conifer dominance likelihood
pred = predict(pcon,newdata=seed.df,type="response",allow.new.levels=TRUE)
pred[pred<0.55]=0
pred[pred!=0]=1

Metrics::accuracy(seed.df$ConRegen,pred)#0.86 #Full model accuracy with 0.55 threshold
Metrics::accuracy(seed.df$ConRegen,rep(1,171))#0.75 #Naive baseline if we always predict conifer

#Now with continuous probabilities
pred = predict(pcon,newdata=seed.df,type="response",allow.new.levels=TRUE)
Metrics::auc(seed.df$ConRegen,pred)#0.85 #Actual AUC. Naive baseline is 0.5 by default

Metrics::logLoss(seed.df$ConRegen,pred)#0.43 #Actual log loss
Metrics::logLoss(seed.df$ConRegen,rep(0.99,171))#1.17 #Naive baseline if we always predict conifer

#Checking spatial autocorrelation in the predictor data
seed_sf <- st_as_sf(seed.df, coords = c("longitude","latitude"))
st_crs(seed_sf) <- st_crs(4326)
seed_sf = st_transform(seed_sf, 3581)#shifting to a projected CRS Canada Lambert

#Read in and mosic the predictor variables
#Take random sample of them for assessing spatial autocorrelation
#This is slow to run
#Needs to be run before the section below reading in the sample
# setwd(#Your path here
#   "/Data/Mod_Pred_Rasts/")
# climlist = list.files(pattern="Seed_Predict")
# climrc = sprc(climlist)
# clim = mosaic(climrc)
# gc()

#This is very slow. Be cautious
#Basically, because there are so many NAs this sample has to be massive
# set.seed(123)
# samp2 = spatSample(clim,
#                    method="regular",
#                    replace = F,
#                    size=10000000,
#                    na.rm=T,
#                    as.points=T)
# writeVector(samp2,
#Your path here
#"/Data/Mod_Pred_Rasts/ClimDat_Sample.shp")


#NOTE YOU CANNOT READ THIS IN IF YOU DON'T RUN THE COMMENTED CODE ABOVE
samp2 = vect(#Your path here
  "/Data/Mod_Pred_Rasts/ClimDat_Sample.shp")

sac1 <- cv_spatial_autocor(x = samp2,
                           column = c("CMI","GDD","Distance_R"),
                           plot=T)
sac1$range#188773.1,
summary(sac1)#Range of autcorr in the predictor vars 188773.1, #590338.6,170351.6

#Create spatial blocks for leave-one spatial block out CV
#Block size is determined by the distance of autcorrelation identified above
set.seed(123)
sb1 <- cv_spatial(
  x = seed_sf,
  column = "ConRegen",#Checks to make sure there are no folds with exclusively only either conifer-dominant or broadleaf-dominant stands
  k = 10, # number of folds
  hexagon=T,
  size = 188773.1, # size of the blocks in meters
  selection = "random", # random blocks-to-fold
  iteration = 50, # find evenly dispersed folds
  progress = FALSE,
  biomod2 = TRUE
)
summary(sb1)

#Here's the blocks
cv_plot(
  cv = sb1, # cv object
  x = seed_sf, # species spatial data
  num_plots = c(1:10) # Plot the folds
)

#Assign folds to the field data
seed.df$folds = sb1$folds_ids

#10-folds spatially blocked
#Cross-validation model performance metrics of tree regen
mae.out = c()
rmsle.out = c()
male.out = c()

foldsu = unique(seed.df$folds)
for(k in 1:length(foldsu)){
  i = foldsu[k]
  train = seed.df[seed.df$folds!=i,]
  test = seed.df[seed.df$folds==i,]
  #REVISE MODEL EACH TIME FOR NEW STATS
  mod = glmmTMB(Total_seed_ha~CMI+GDD+Distance_Reburnedge_m+years_post+(1|project_id),
                ziformula=~Distance_Reburnedge_m+GDD,
                data=train,
                family=poisson())
  pred = predict(mod,newdata=test,type="response",allow.new.levels=TRUE)#re.form=NA
  rmsle.1 = rmsle(test$Total_seed_ha,pred)
  mae.1 = mae(test$Total_seed_ha,pred)
  male.1 = mean(abs(log10((test$Total_seed_ha+1)/pred+1)))
  
  mae.out[k] = mae.1
  rmsle.out[k] = rmsle.1
  male.out[k] = male.1
}
mean(mae.out)#28064.12
mean(rmsle.out)#2.46
mean(male.out)#0.24

#Now for the %conifer model. We will assess this one with balanced accuracy of our threshold value (0.55)
library(cutpointr)
seed.df$Con_Pred = predict(pcon,type="response",newdata = seed.df)
point=cutpointr(seed.df,Con_Pred,ConRegen,tol_metric = 0.05,method = maximize_metric, 
                metric = youden,break_ties = min,boot_runs=1000,use.midpoints=T)
t.thresh = point$optimal_cutpoint #Optimal cutpoint = 0.55
t.thresh
plot(point)

set.seed(123)
foldsu = unique(seed.df$folds)
acc.out = c()
auc.out = c()
ll.out = c()
for(k in 1:length(foldsu)){
  i = foldsu[k]
  train = seed.df[seed.df$folds!=i,]
  set.seed(123)
  seed.over = ovun.sample(ConRegen~., data = train, method = "both", N = round(nrow(train)*1.169591))
  seed.over = seed.over$data
  test = seed.df[seed.df$folds==i,]
  #fitting model with new randomly up/down sampled dataset from remaining folds
  #Change model to compare
  mod = glmmTMB(ConRegen~CMI+Distance_Reburnedge_m,#+(1|project_id),
                data=seed.over,
                family=binomial)
  #Test only on holdout fold
  pred = predict(mod,newdata=test,type="response",allow.new.levels=TRUE)
  pred[pred<0.55]=0
  pred[pred!=0]=1
  accuracy.1 = Metrics::accuracy(test$ConRegen,pred)
  pred = predict(mod,newdata=test,type="response",allow.new.levels=TRUE)
  ll.1 = Metrics::logLoss(test$ConRegen,pred)
  auc.1 = Metrics::auc(test$ConRegen,pred)
  acc.out[k] = accuracy.1
  auc.out[k] = auc.1
  ll.out[k] = ll.1
}
mean(acc.out)#0.83
mean(auc.out)#0.8034607
mean(ll.out)#0.4463502
