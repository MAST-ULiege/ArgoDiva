source("../../Utils/ArgoSelect.R")

dataDir   = "../../Datas/Coriolis_PAR/"
selectCriterium0<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"*\\.nc"),#[1:5],
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("DOWNWELLING_PAR"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
fulldf <- ArgoSelect(selectCriterium0)

fulldf$variable[which(fulldf$variable=="DOWNWELLING_PAR")]<-"PAR"




