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

# To now about defined function you could use , look at :
source("ArgoVertFunctions.R")

# To add new function, use :  (using same structure as in the previous file)
source("ArgoVertFunctions_USER.R")

# list of the diagnostic to extract from Argos
flist <- list("MaxTem","MaxSal")

fdf   <- ArgoExtract(fulldf, flist)

##################################
# Paremeters for diva input file #
##################################



