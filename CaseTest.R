# Case Test example for the use of the ArgoDiva toolbox
# A. Capet 02/02/2018

# Load package and function sources
source("ArgoLoad.R")

################################
# Criterium for Argo selection #
################################

selectCriterium<-list(
  # List of individual Argo netcdf files. 
  # *TODO* provide an option to give only the name of a repertory
  dataDir   = "~/Desktop/DOCS/TEACHING/Florient/ArgoData/",
  filenames = c("6901866_Mprof.nc"), 
  
  # List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL"), 
  
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)

fulldf <- ArgoSelect(selectCriterium)

########################
# Function definitions #
########################

source("ArgoVertFunctions.R")






##################################
# Paremeters for diva input file #
##################################



