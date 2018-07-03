# Case Test example for the use of the ArgoDiva toolbox
# A. Capet 02/02/2018

# Load package and function sources
source("ArgoLoad.R")

################################
# Criterium for Argo selection #
################################

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/PUBLISHED/VOX/ARGODATACIL2/"

selectCriterium<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  # filenames = c("6901866_Mprof.nc"),
  # filenames = c("argo-profiles-5902291.nc","argo-profiles-6901960.nc"),
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"*.nc")[1:5],
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
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
flist <- list("DepthMinTem")

fdf   <- ArgoExtract(fulldf, flist)

###########
# Display #
###########

#Visualization of distribution density for different period (showing different approach)
ArgoDisplay(fdf,selectCriterium,flist)

##################################
# Paremeters for diva input file #
##################################



