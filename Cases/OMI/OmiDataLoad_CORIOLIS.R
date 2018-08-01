source("../../Utils/ArgoSelect.R")

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2010-2014/"
selectCriterium0<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
#fulldf0 <- ArgoSelect_CMEMS(selectCriterium0, datasource = "Coriolis")
fulldf0 <- ArgoSelect(selectCriterium0)

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2014-2016/"
selectCriterium1<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
#fulldf1 <- ArgoSelect_CMEMS(selectCriterium1, datasource = "Coriolis")
fulldf1 <- ArgoSelect(selectCriterium1)

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2017/"
selectCriterium2<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  # filenames = c("6901866_Mprof.nc"),
  # filenames = c("argo-profiles-5902291.nc","argo-profiles-6901960.nc"),
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  #filenames = Sys.glob(paste0(dataDir,"argo-profiles-*.nc")),
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
fulldf2 <- ArgoSelect(selectCriterium2)

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2018/"
selectCriterium3<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  # filenames = c("6901866_Mprof.nc"),
  # filenames = c("argo-profiles-5902291.nc","argo-profiles-6901960.nc"),
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  #filenames = Sys.glob(paste0(dataDir,"argo-profiles-*.nc")),
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
fulldf3 <- ArgoSelect(selectCriterium3)

fuldf<-rbind(fulldf0,fulldf1,fulldf2,fulldf3)


