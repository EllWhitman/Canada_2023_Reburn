## **2023 wildfires expose Canadian forests to transformative state change**


## *Overview*
The Canada_2023_Reburn repo replicates the analysis presented in the manuscript ‘2023 wildfires expose Canadian forests to transformative state change’. The code is written to be run in R and requires supporting packages which are referenced at the start of each script.
The scripts generate maps of overlaps between mapped polygon-based fire perimeters and analyse the size, frequency, and total area reburned over time. The user can set the window of time between fires to define a reburn (e.g., 10, 20, 30 years between fires). The scripts then identify core area (far from seed sources) of reburns, and predict likely post-fire vegetation community trajectories within reburned areas, as a function of climate (aridity and temperature) and distance to reburn edge as a proxy for tree seed sources. The code also produces figures displaying the analysis, and summary statistics. For assistance with this repository or scripts please contact Ellen Whitman (ellen [dot] whitman [at] nrcan-rncan.gc.ca).
**This analysis is presented in a manuscript currently in review.**

## *Repo Contents*
[R](https://github.com/EllWhitman/Canada_2023_Reburn/tree/main/R): r code

## *System requirements*
Hardware requirements:
A computer with adequate RAM to manipulate spatial data is necessary for analysing large areas, such as the entire boreal zone of Canada. The runtimes reported are from a computer with a 8 cores with 2.4 GHz per core, and 64 GB of RAM. 

## *Software requirements*
Scripts were written and run on MacOS Sequoia 15.6.1 in RStudio, using R version 4.5.2. Users require R to be installed on their computer in order to run the scripts. https://cran.r-project.org/  

## *How to use*
Download the scripts, after having installed R, and open them in R or an R development environment where they can be run. Wherever there is a comment in the code saying ‘Your path here’ edit the code to direct to the user’s own local path to locate the file. The scripts require dependencies of certain packages, which are called at the start of each script. You may need to install these packages if you have not already. The scripts are to be followed in numbered order, beginning with 1.

## *Run time*
The first two scripts, which manipulate national spatial datasets in a loop, may require several hours (e.g., an entire work day) to run depending on the area of interest and duration of the time series selected by the user. I recommend selecting a test time period (e.g., 3 consecutive years) to trial the code, determine if the selected study area is reasonable for the computing power available, and to gauge what the runtime of a longer time series will be. The later scripts (numbers 3, 5 - 8) each run in approximately 1 minute to 20 minutes. The exception to this runtime is script 4 where spatial predictions of post-fire recruitment density and likelihood of conifer dominance are made from models, which also requires several hours depending on the extent of the area of interest. Because the field data used to fit these models comes from a variety of sources, not all of which are public, script 4 is bypassed by downloading the test data from the data folder. This data consists of spatial predictions from the two models fitted in script 4 within a small test area in northwestern Canada.
