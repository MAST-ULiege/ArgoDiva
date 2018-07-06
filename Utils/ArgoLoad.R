# libraries required for the ArgoDiva Toolbox 
library(ncdf4)
library(ggplot2)
library(plyr)
library(reshape2)
library(chron)

basedir <- '/home/arthur/Desktop/DOCS/TEACHING/Florient/ArgoDiva'


source(paste0(basedir,'/',"Utils/ArgoSelect.R"))
source(paste0(basedir,'/',"Utils/ArgoExtract.R"))
source(paste0(basedir,'/',"Utils/ArgoDisplay.R"))
source(paste0(basedir,'/',"Utils/ArgoPrepareforDiva.R"))

