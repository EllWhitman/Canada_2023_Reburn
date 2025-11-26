#######################
##Mapping Area Reburned in the Canadian Boreal Forest in 2023, and prior years
#Analysis of climate and core area within 2023 short-interval reburned areas
#Written by Ellen Whitman
#ellen.whitman@nrcan-rncan.gc.ca
#Last Edit Nov 25, 2025
#Wherever 'Your path here' appears as a comment, please replace with the correct pathway to your local file

#Read in predicted trajectory within reburns data. Each row = 1ha
climdat = read.csv(#Your path here
  "/Data/Regen_Trajectories_Data.csv")

table(climdat$Trajectory)#Total area of each trajectory, as originally modelled

#Create new trajectories, where definition of sparse is more strict
climdat$Trajectory100 = NA
climdat[climdat$Total_Stems_ha<=100,"Trajectory100"] = "Sparse"#100 stems/ha = sparse
climdat[climdat$Total_Stems_ha>100&climdat$PConifer<0.56,"Trajectory100"]="Broadleaf"
climdat[is.na(climdat$Trajectory100),"Trajectory100"]="Conifer"

table(climdat$Trajectory100)#Total area of each trajectory, when sparse is 100 stems/ha

#Create new trajectories, where definition of sparse is less strict (corresponds to stem density of mature stand)
climdat$Trajectory5000 = NA
climdat[climdat$Total_Stems_ha<=5000,"Trajectory5000"] = "Sparse"#5000 stems/ha = sparse
climdat[climdat$Total_Stems_ha>5000&climdat$PConifer<0.56,"Trajectory5000"]="Broadleaf"
climdat[is.na(climdat$Trajectory5000),"Trajectory5000"]="Conifer"

table(climdat$Trajectory5000)#Total area of each trajectory, when sparse is 5000 stems/ha
