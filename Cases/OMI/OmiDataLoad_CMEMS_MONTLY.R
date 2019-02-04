dataDirBase <- paste0(basedir,'/Datas/CMEMS_MONTHLY/')
dataDirList <- list.dirs(dataDirBase)

fuldflist <- list()

for (i in seq_along(dataDirList)) {
selectCriterium0<-list(
    dataDir   = dataDirList[[i]],
    ## List of individual Argo netcdf files. 
    ########################################
    ## To get all files in the directory use 
    filenames = list.files(dataDir,"GL_PR_PF_.*\\.nc"),#[1:5],
    ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
    varList      = c("TEMP","PSAL","DOX2"), 
    # This flag states wether we retain only records for which all variables are presents
    presentInAll = TRUE
  )


fuldf <- ArgoSelect_CMEMS(selectCriterium0, datasource = "CMEMS")
fuldf$variable[which(fuldf$variable=="DOX2")]<-"DOXY"
fuldflist[[i]] <- fuldf
}

fuldf <- do.call(fuldflist, rbind)