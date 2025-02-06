## **Extreme 2023 wildfires driving a biome shift in the Canadian boreal forest**


## *Overview*
The Canada_2023_Reburn repo replicates the analysis presented in the manuscript ‘Extreme 2023 wildfires driving a biome shift in the Canadian boreal forest’. The code is written to be run in R and requires supporting packages which are referenced at the start of each script.
The scripts generate maps of overlaps between mapped polygon-based fire perimeters and analyse the size, frequency, and total area reburned over time. The user can set the window of time between fires to define a reburn (e.g., 10, 20, 30 years between fires). The scripts then identify core area (far from seed sources) of reburns, and mapped reburns are compared to raster maps of climate normals to identify climate-stressed reburned areas. The code also produces figures displaying the analysis, and summary statistics. For assistance with this repository or scripts please contact Ellen Whitman (ellen [dot] whitman [at] nrcan-rncan.gc.ca).
This analysis is presented in a manuscript currently in review.

## *Repo Contents*
[R](https://github.com/EllWhitman/Canada_2023_Reburn/tree/main/R): r code

## *System requirements*
Hardware requirements
A computer with adequate RAM to manipulate spatial data is necessary for analysing large areas, such as the entire boreal zone of a Canada. The runtimes reported are from a computer with a 8 cores with 2.4 GHz per core, and 64 GB of RAM. 

## *Software requirements*
Scripts were written and run on MacOS Sonoma in RStudio, using R version 4.4.2. Users require R to be installed on their computer in order to run the scripts. https://cran.r-project.org/  

## *How to use*
Download the scripts, after having installed R, and open them in R or an R development environment where they can be run. Wherever there is a comment in the code saying ‘Your path here’ edit the code to direct to the user’s own local path to locate the file. The scripts require dependencies of certain packages, which are called at the start of each script. You may need to install these packages, if you have not already. The scripts are to be followed in numbered order, beginning with 1.

## *Runtime*
The first two scripts, which manipulate national spatial datasets in a loop, may require several hours (e.g., an entire work day) to run, depending on the area of interest and duration of the time series selected by the user. I recommend selecting a test time period (e.g., 3 consecutive years) to trial the code, determine if the selected study area is reasonable for the computing power available, and to gauge what the runtime of a longer time series will be. The later scripts (numbers 3 – 5) each run in approximately 1 minute, 20 minutes, and <1 minute, respectively.
