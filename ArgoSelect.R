# This function extract Argo profiles based on user input criterium
# 
# Arguments (in) :
# 
# Value (out) :
#  * A dataframe with columns :
#         - value 
#         - qc
#         - alevel   : the index of the Argo Level
#         - depth    : the depth [m]
#         - aprofile : the index of the profile 
#         - variable : the name of the variable
#         - juld     : The julian day 
#         - lon      : Longitude [decimal °]
#         - lat      : Latitude [decimal °]

ArgoSelect <- function(selectCriterium){
  with(selectCriterium,{
    
    argodflist<-lapply(as.list(filenames),function(file){
      
      #Opening the file in a open-only mode
      ncfile   <<- nc_open(paste0(dataDir,file), write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
      
      #Dimensions
      N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
      N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
      juld     <- ncvar_get(ncfile,"JULD")
      pres     <- ncvar_get(ncfile,"PRES")
      lon      <- ncvar_get(ncfile,"LONGITUDE")
      lat      <- ncvar_get(ncfile,"LATITUDE")
      
      FloatInfo<-list(N_PROF=N_PROF,
                      N_LEVELS=N_LEVELS,
                      juld = juld,
                      pres=pres,
                      lon=lon,
                      lat=lat)
      
      dflist <-lapply (varList, function (V){
        vstringadj<-paste0(V,'_ADJUSTED')
        vardf <- ExtractVar(vstringadj,FloatInfo)
        if (all(is.na(vardf)) == TRUE) {   # is "= TRUE" needed ? 
          vardf <- ExtractVar(V,FloatInfo)
        }  
      })
      
      dftot <- do.call(rbind, dflist)

      id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
      
      dftot$Platform <- unique(id)
      dftot$day      <- month.day.year(dftot$juld,c(1,1,1950))$day
      dftot$month    <- month.day.year(dftot$juld,c(1,1,1950))$month
      dftot$year     <- month.day.year(dftot$juld,c(1,1,1950))$year
      
      return(dftot)
     })
    return (do.call(rbind,argodflist))
  })
}

# This function should return a dataframe for variable with value, qc, iprofile and ilevel
ExtractVar<-function(Var,FloatInfo){
  with(FloatInfo,{
    lvar             <- ncvar_get(ncfile,Var)
    lvar_qc          <- ncvar_get(ncfile,paste0(Var,"_QC"))
    lvar_qctab<-llply(lvar_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab<-do.call(cbind,lvar_qctab)
    
    # making dataframes, removing the NANs  
    alevels<-1:N_LEVELS
    d<-ldply(as.list(1:N_PROF),function(iprof){
      indexes<- !(is.na(lvar[,iprof])|is.na(lvar[,iprof]))
      if(sum(indexes)==0){
        return (data.frame())
      }
      
      data.frame(value = lvar[indexes,iprof],
                 qc      = as.integer(lvar_qctab[indexes,iprof]),
                 alevel  = alevels[indexes],
                 depth   = pres[indexes,iprof],
                 aprofile = iprof,
                 variable = Var)
    })
    
    d$juld <-juld[d$aprofile]
    d$lon  <-lon[d$aprofile]
    d$lat  <-lat[d$aprofile]
    
    return(d=d)
  })
}  

