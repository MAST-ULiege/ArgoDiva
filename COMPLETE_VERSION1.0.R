library(ncdf4)
library(ggplot2)
library(plyr)
library(reshape2)
library(chron)
library(viridis)
library(TTR)
library(minpack.lm)
library(nls2)
library(nlstools)
library(data.table)
library(gsw)
library(gridExtra)
library(car)
library(oce)

#/!!!!!!!!!!!!!!!!!!!!!!!\     CHECK BEFORE USE     /!!!!!!!!!!!!!!!!!!!!\

check_data <- "FLUO"
check_mode <- "ADJUSTED"
check_sigma <- "TEOS10"
check_mld <- "BIRA"

#FLUORESCENCE PARAMETERS AVAILABLE IN METADATA FILES FROM EACH FLOATS
param_6900807 <- c(0.0072, 47) #SCALE_CHLA, DARK_CHLA
param_6901866 <- c(0.0073, 49)
param_7900591 <- c(0.0121, 48)
param_7900592 <- c(0.0121, 50)
fluo_param <- list(param_6900807, param_6901866, param_7900591,
                   param_7900592)

#CHLA-equipped ARGO floats
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

# FUNCTION TO BE USED
ExtractVar <- function(Var,FloatInfo){
  with(FloatInfo,{
    # This function should return a dataframe for variable with value, qc, iprofile and ilevel
    lvar             <- ncvar_get(ncfile,Var)
    lvar_qc          <- ncvar_get(ncfile,paste0(Var,"_QC"))
    lvar_qctab       <- llply(lvar_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab<-do.call(cbind,lvar_qctab)
    
    # making dataframes, removing the NANs  
    alevels <- 1:N_LEVELS
    d <- ldply(as.list(1:N_PROF),function(iprof){
      indexes <- !(is.na(lvar[,iprof])|is.na(lvar[,iprof]))
      if(sum(indexes) == 0){
        return (data.frame())
      }
      
      data.frame(value    = lvar[indexes,iprof],
                 qc       = as.integer(lvar_qctab[indexes,iprof]),
                 alevel   = alevels[indexes],
                 depth    = pres[indexes,iprof],
                 aprofile = iprof,
                 variable = Var)
    })
    
    d$juld <- juld[d$aprofile]
    d$lon  <- lon[d$aprofile]
    d$lat  <- lat[d$aprofile]
    
    return(d=d)
  })
}  

#MEDIAN FILTER OVER 5 POINTS
mmed <- function(x,n=5){runmed(x,n)}

#NORMALIZATION (FOR FLUO DATA BEFORE GAUSSIAN FITTING FOR INSTANCE)
normalize <- function(data){
  data <- (data - min(data, na.rm = T))/(max(data, na.rm = T) - min(data, na.rm = T))
}

#Get each profile from ARGO data -----
profiles <- ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
  #Dimensions
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  
  FloatInfo <- list(N_PROF=N_PROF,
                    N_LEVELS=N_LEVELS,
                    juld=juld,
                    pres=pres,
                    lon=lon,
                    lat=lat)
  
  if(check_mode == "ADJUSTED"){
  #FOR NPQ AUTOMATIC CORRECTION (SEE SCHMECHTIG 2014)
  #CHLA_ADJUSTED IS A RT CORRECTION (RT : REAL TIME)
    chladf <- ExtractVar("CHLA_ADJUSTED", FloatInfo) 
  }else{
    chladf <- ExtractVar("CHLA", FloatInfo)
  }
  # CHLA = (FLUO_CHLA - DARK_CHLA) * SCALE_CHLA
  # => FLUO_CHLA = (CHLA / SCALE_CHLA) + DARK_CHLA
  
  #fluorescence conversion, respectively to each float
  if (check_data == "FLUO"){
    if (file == "6900807_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[1]][1]) + fluo_param[[1]][2] 
    }
    else if (file == "6901866_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[2]][1]) + fluo_param[[2]][2] 
    }
    else if (file == "7900591_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[3]][1]) + fluo_param[[3]][2] 
    }
    else if (file == "7900592_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[4]][1]) + fluo_param[[4]][2] 
  }}
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  data.frame(depth    = chladf$depth,
             juld     = chladf$juld,
             fluo     = chladf$value,
             qc       = chladf$qc,
             day      = month.day.year(chladf$juld,c(1,1,1950))$day,
             month    = month.day.year(chladf$juld,c(1,1,1950))$month,
             year     = month.day.year(chladf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = chladf$lon,
             lat      = chladf$lat,
             Platform = as.numeric(unique(id)),
             type     = "Argo")
})

#REMOVE BAD DATA (I.E. QC = 4 FOR NEGATIVE SPIKES, JUMPS, ETC.)
#DELETION OF DEPTH WITH BAD DATA
profiles <- profiles[-(which(profiles$qc == 4)),] 

#CREATION OF PROFILE IDs & REORDER DATA FRAME ACCORDING TO IT
profiles <- transform(profiles,id=as.numeric(factor(juld)))
profiles <- profiles[order(profiles$id),]

#MEDIAN FILTER TO REMOVE POSITIVE SPIKES (THEY MAY RETAIN INFORMATION BUT STRONG PIKES 
#WILL NOT BE RETAINED FOR THIS STUDY)
smoothed_fluo <- ldply(as.list(unique(profiles$id)), function(i){
  tmp <- profiles[profiles$id == i,]
  tmp <- mmed(tmp$fluo, 5)
  data.frame(smoothed_fluo = tmp)
})

#REPLACE FLUO WITH SMOOTHED FLUO
profiles[,3] <- smoothed_fluo

#Construction of an ARGO-only dataframe -> first guess value for Zmax [NLS fit]
argodf <- ldply(as.list(filename), function(file){
  
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
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
  
  if(check_mode == "ADJUSTED"){
    chladf <- ExtractVar("CHLA_ADJUSTED", FloatInfo) 
  }else{
    chladf <- ExtractVar("CHLA", FloatInfo)
  }
  
  if (check_data == "FLUO"){
    if (file == "6900807_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[1]][1]) + fluo_param[[1]][2] 
    }
    else if (file == "6901866_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[2]][1]) + fluo_param[[2]][2] 
    }
    else if (file == "7900591_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[3]][1]) + fluo_param[[3]][2] 
    }
    else if (file == "7900592_Mprof.nc"){
      chladf$value <- (chladf$value / fluo_param[[4]][1]) + fluo_param[[4]][2] 
    }}
  
  subchladf <- subset(chladf,select=c("depth","juld","value","qc","lon","lat"))
  colnames(subchladf)[which(colnames(subchladf)=="value")]<-"fluo"
  
  chladf <- ddply(subchladf,~juld,summarize,
                  qc = qc[which.max(fluo)],
                  depthmax = depth[which.max(fluo)],
                  maxvalue = fluo[which.max(fluo)],
                  integration = sum(fluo),
                  bottomdepth = max(depth),
                  lon=mean(lon),
                  lat=mean(lat))
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  data.frame(depthmax        = chladf$depthmax,
             bottom          = chladf$bottomdepth,
             qc              = chladf$qc,
             maxvalue        = chladf$maxvalue,
             integratedvalue = chladf$integration,
             juld            = chladf$juld,
             day             = month.day.year(chladf$juld,c(1,1,1950))$day,
             month           = month.day.year(chladf$juld,c(1,1,1950))$month,
             year            = month.day.year(chladf$juld,c(1,1,1950))$year,
             DOY             = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), 
                                                   format ="%j")),#Day Of Year
             lon             = chladf$lon,
             lat             = chladf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
})

#IDs
argodf <- transform(argodf,id=as.numeric(factor(juld)))
argodf <- argodf[order(argodf$juld),]

#DENSITY PROFILES
density_profiles <- ldply(as.list(filename),function(file){
  
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
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
  
  #NO ADJUSTED VALUE FOR PSAL AND TEMP FOR EACH FLOAT
  tempdf <- ExtractVar("TEMP",FloatInfo)
  tempdf$value[which(tempdf$qc == 4 | tempdf$qc == 3)] <- NA
  psaldf <- ExtractVar("PSAL",FloatInfo)
  psaldf$value[which(psaldf$qc == 4 | psaldf$qc == 3)] <- NA
  
  if(check_sigma =="TEOS10"){
  #DENSITY ANOMALY BASED ON TEOS-10
  psal <- gsw_SA_from_SP(psaldf$value,psaldf$depth,psaldf$lon,psaldf$lat)
  temp <- gsw_CT_from_t(psal,tempdf$value,tempdf$depth)
  sigma <- gsw_sigma0(psal,temp)
  }else{
    sigma <- swSigmaTheta(psaldf$value,tempdf$value,tempdf$depth)
  }
  
  data.frame(sigma = sigma,
             depth = tempdf$depth, 
             juld  = tempdf$juld, 
             lon   = tempdf$lon,
             lat   = tempdf$lon,
             day   = month.day.year(tempdf$juld,c(1,1,1950))$day,
             month = month.day.year(tempdf$juld,c(1,1,1950))$month,
             year  = month.day.year(tempdf$juld,c(1,1,1950))$year,
             DOY   = as.integer(strftime(as.Date(tempdf$juld,origin = '1950-01-01'), 
                                                 format ="%j"))#Day Of Year
  )
})

#IDs
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))
density_profiles <- density_profiles[order(density_profiles$id),]

#Get profiles id that do not answer to some criteria
crit1 <- argodf[(argodf$bottom < 100),]#TOO SHALLOW PROFILES
crit2 <- argodf[(argodf$depthmax >= 100),]#NOT USUAL AT ALL ABOVE 100M
crit3 <- argodf[(argodf$depthmax == argodf$bottom),]#STUCK PROFILES AT SAME DEPTH
criteriadf <- rbind(crit1, crit2, crit3)
criteriadf <- unique(criteriadf)

#remove profiles that have been categorized as bad profiles
for (i in 1:length(criteriadf$id)){
  argodf <- argodf[!(argodf$id == criteriadf$id[i]),]
  profiles <- profiles[!(profiles$id == criteriadf$id[i]),]
  density_profiles <- density_profiles[!(density_profiles$id == criteriadf$id[i]),]
}

#REMOVE 2013 ET 2018 (because not FULL year)
argodf <- argodf[!(argodf$year == 2013 | argodf$year == 2018),]
profiles <- profiles[!(profiles$year == 2013 | profiles$year == 2018),]
density_profiles <- density_profiles[!(density_profiles$year == 2013 | density_profiles$year == 2018),]

#new id before gaussian elimination
profiles <- transform(profiles,id=as.numeric(factor(juld)))
argodf <- transform(argodf,id=as.numeric(factor(juld)))
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))

#REMOVE DEPTH WITH NA VALUES IN DENSITY PROFILES (WE MAY LOSE FULL PROFILES WITH THAT)
#density_profiles <- density_profiles[-(which(is.na(density_profiles$sigma) == TRUE)),]

#COMPUTE MLD
if(check_mld == "BIRA"){
  sigma_criteria <- 0.125
  depth_ref <- 3
}else if(check_mld == "ARGO"){
  sigma_criteria <- 0.03
  depth_ref <- 10
}else{#PERSONAL CRITERIA BASED ON ARTHUR'S ADVICE AND USED BY XING (but assume a min MLD of 15 m...)
  sigma_criteria <- 0.03
  depth_ref <- 15
}


MLDdf <- ldply(as.list(1:length(unique(density_profiles$id))), function(i){
  
  tmp <- density_profiles[density_profiles$id == i,]
  
  if(all(is.na(tmp$sigma)) == TRUE){
    MLD <- NA
    sigma_surface <- NA
  }else if(length(tmp$sigma[!is.na(tmp$sigma)==TRUE])>=2){
    sigma_surface<-NA
    sigma_surface <- approx(tmp$depth,tmp$sigma,depth_ref)$y
  }
  
    if(is.na(sigma_surface) == FALSE){
      MLD <- max(tmp$depth[tmp$sigma <= (sigma_surface + sigma_criteria)], na.rm = T)
    }else{
      MLD <- NA
    }
  
  data.frame(sigma_surface = sigma_surface, MLD = MLD, juld = tmp$juld[1], day = month.day.year(tmp$juld[1],c(1,1,1950))$day,
             month = month.day.year(tmp$juld[1],c(1,1,1950))$month,
             year = month.day.year(tmp$juld[1],c(1,1,1950))$year,
             DOY = as.integer(strftime(as.Date(tmp$juld[1],origin = '1950-01-01'), 
                                       format ="%j")))#Day Of Year
  
})

MLDdf <- transform(MLDdf,id=as.numeric(factor(juld)))

quenching_correction <- function(fluo,depth,MLD) {
  if(is.na(MLD) == FALSE){
    f <- fluo[!is.na(fluo) & depth <= MLD]
    d <- depth[!is.na(fluo) & depth <= MLD]
    zMax <- d[which.max(f)]
    Max <- max(f)
    Corfluo <- fluo
    #Criteria from Schmechtig et al. 2014
    if(!is.na(MLD) & min(f[d<=zMax])<=(0.9*Max)) Corfluo[depth<=zMax] <- Max
    return(Corfluo)
  }else{
    return(fluo)
  }
}

if(check_mode != "ADJUSTED"){
  
  #NEED TO FIND A WAY TO DETERMINE ALL MLD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  remove_uncertain_profiles <- which(is.na(MLDdf$rho_surf))
  for (i in 1:length(remove_uncertain_profiles)){
    profiledf <- profiledf[!(profiledf$id == remove_uncertain_profiles[i]),]
    MLDdf <- MLDdf[!(MLDdf$id == remove_uncertain_profiles[i]),]
  }
  profiledf <- transform(profiledf,id=as.numeric(factor(juld)))
  MLDdf <- transform(MLDdf,id=as.numeric(factor(juld)))
  
  #NOTE que si la première profondeur est fort profonde (du coup) -> Quenching est négligeable 
  # et le fait que pas de correction (si MLD = NA) ne crée pas de 'FAUX' DCM -> à regarder
  
  ##################################################################################################
  
  fluo_NPQ <- ldply(as.list(1:length(unique(profiles$id))), function(i){
    tmp <- profiles[profiles$id == i,]
    correction <- quenching_correction(tmp$fluo, tmp$depth, MLDdf$MLD[i])
    data.frame(fluo_NPQ = correction)
})
  profiles$fluo <- fluo_NPQ$fluo_NPQ
}


#Following Mignot et al. 2011
fgauss <- function(x, Fsurf, Zdemi, Fmax, Zmax, dz){
  Fsurf*exp((-log(2)/Zdemi)*x) + Fmax*exp(-(x-Zmax)^2/dz^2)
}

fsigmoid <- function(x, Fsurf, Zdemi, s){
  Fsurf*(1/(1+exp((Zdemi-x)*s)))
}

gaussiandf <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- profiles[profiles$id==i,]
  depthindex <- which.min(tmp$depth <= 100)
  
  if(check_data == "FLUO"){
    tmp$fluo <- normalize(tmp$fluo)
  }
  
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
  maxindex <- which.max(tmp$fluo)

  #Parameters estimation
  nonNAindex <- which(!is.na(tmp$fluo))
  firstnonNA <- nonNAindex[1]
  Fsurf <- tmp$fluo[firstnonNA]
  Zmax <- tmp$depth[maxindex]
  if(Zmax > 100){
    Zmax <- 50
  }
  Fmax <- tmp$fluo[maxindex]
  Zdemi <- tmp$depth[which(tmp$fluo <= Fsurf/2)[2]]
  if(is.na(Zdemi)){#code à la bourrain
    indexZ <- which(tmp$fluo < Fsurf)
    Zdemi <- tmp$depth[indexZ[which.max(indexZ > maxindex)]]
  }
  dz <- 15
  
  res <- tryCatch(nlsLM(y ~ fgauss(x, Fsurf, Zdemi, Fmax, Zmax, dz),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                        data=tab, control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  
  if(inherits(res, "error")){
    dz <- 10
    res <- tryCatch(nlsLM(y ~ fgauss(x, Fsurf, Zdemi, Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res, "error")){
    dz <- 5
    res <- tryCatch(nlsLM(y ~ fgauss(x, Fsurf, Zdemi, Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res, "error")){
    dz <- 20
    res <- tryCatch(nlsLM(y ~ fgauss(x, Fsurf, Zdemi, Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  
  
  data.frame(Fsurf = coef(res)["Fsurf"], Zdemi = coef(res)["Zdemi"],
             Fmax = coef(res)["Fmax"], Zmax = coef(res)["Zmax"], dz = coef(res)["dz"], id=i,
             juld=tmp$juld[1], file=tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
  
})

sigmoidf <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- profiles[profiles$id==i,]
  depthindex <- which.min(tmp$depth <= 100)
  
  if(check_data == "FLUO"){
    tmp$fluo <- normalize(tmp$fluo)
  }
  
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
  maxindex <- which.max(tmp$fluo)
  
  #Parameters estimation
  nonNAindex <- which(!is.na(tmp$fluo))
  firstnonNA <- nonNAindex[1]
  Fsurf <- tmp$fluo[firstnonNA]
  Zdemi <- tmp$depth[which(tmp$fluo <= Fsurf/2)[1]]
  if(is.na(Zdemi)){
    indexZ <- which(tmp$fluo < Fsurf)
    Zdemi <- tmp$depth[indexZ[which.max(indexZ > maxindex)]]
  }
  i1 <- which.min(tmp$depth <= Zdemi - 5)
  i2 <- which.min(tmp$depth <= Zdemi + 5)
  d1 <- tmp$depth[i1]
  d2 <- tmp$depth[i2]
  s <- -(Fsurf/2)/(d2-d1)
  
  res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                        data=tab, control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  if(inherits(res, "error")){
    s <- -0.1
    res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res, "error")){
    s <- -0.01
    res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res, "error")){
    s <- -1
    res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  data.frame(Fsurf = coef(res)["Fsurf"], Zdemi = coef(res)["Zdemi"], s=coef(res)["s"], id=i,
             juld=tmp$juld[1], file=tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
  
})


# i <- 1
# tmp <- profiles[profiles$id==i,]#gappy data, no interpolation
# if(check_data == "FLUO"){
#   tmp$fluo <- normalize(tmp$fluo)
# }
# depthindex <- which.min(tmp$depth <= 100)
# tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
# x <- tab$x
# plot(y~x, data=tab, type="l", lwd = 2)
# lines(x = tab$x, fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
#                           sigmoidf$s[i]),col=2, add=T, xlim=range(tab$x), lwd = 2)
# lines(x = tab$x, fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
#                         gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i]),
#       col=4, add=T, xlim=range(tab$x), lwd = 2)
# i <- i + 1

Rcoefdf <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  tmp <- profiles[profiles$id==i,]
  depthindex <- which.min(tmp$depth <= 100)
  
  if(check_data == "FLUO"){
    tmp$fluo <- normalize(tmp$fluo)
  }
  
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
  x <- tab$x
  sigmoid <- fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
                      sigmoidf$s[i])
  
  gauss <- fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
                  gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i])
  
  mean_data <- mean(tab$y, na.rm=T)
  ss_tot <- sum((tab$y - mean_data)^2, na.rm=T)
  ss_res_sigmoid <- sum((sigmoid - tab$y)^2, na.rm = T)
  ss_res_gaussian <- sum((gauss - tab$y)^2, na.rm = T)
  rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)
  rcoef_gaussian <- 1-(ss_res_gaussian/ss_tot)
  data.frame(rcoef_sigmoid = rcoef_sigmoid, rcoef_gaussian = rcoef_gaussian, id = tmp$id[1], juld = tmp$juld[1],
             file = tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
})

#PROFILES CLASSIFICATION
R_crit <- 0.8 #80% OF EXPLAINED VARIANCE (90% FOR NAVARRO)

classification <- ldply(as.list(1:length(unique(profiles$id))), function(i){

  tmp <- gaussiandf[gaussiandf$id==i,]
  
  if(Rcoefdf$rcoef_sigmoid[i] & Rcoefdf$rcoef_gaussian[i] < R_crit){
    classif <- "other"
    score <- "NA"#
  } 
  else if(Rcoefdf$rcoef_sigmoid[i] >= Rcoefdf$rcoef_gaussian[i]){
    classif <- "sigmoid"
    score <- Rcoefdf$rcoef_sigmoid[i]
  }
  else{
    classif <- "gaussian"
    score <- Rcoefdf$rcoef_gaussian[i]
  }
  
  data.frame(classif = classif, id = i, juld = tmp$juld[1], file = tmp$file[1],
             depthmax = tmp$Zmax[1], score = score, lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
})

#count profiles
ngauss <- length(which(classification$classif == "gaussian"))/length(unique(profiles$id))*100
nsigmoid <- length(which(classification$classif == "sigmoid"))/length(unique(profiles$id))*100
nother <- length(which(classification$classif == "other"))/length(unique(profiles$id))*100

#gaussian profile (classified) before a new criteria elimination
gaussprofiles <- classification[classification$classif == "gaussian",]

#remove profiles where depthmax (DCM) < 1m [Navarro]
gaussprofiles <- gaussprofiles[gaussprofiles$depthmax > 1,]

#NEW IDs
gaussprofiles <- transform(gaussprofiles,id=as.numeric(factor(juld)))

#GET 'PSEUDO' GAUSSIAN PROFILES
gaussianprofdf <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- profiles[profiles$juld == gaussprofiles$juld[i],]
  data.frame(depth = tmp$depth, juld = tmp$juld, fluo = tmp$fluo, qc = tmp$qc,
             day = tmp$day, month = tmp$month, year = tmp$year, DOY = tmp$DOY,
             Platform = tmp$Platform, lon = tmp$lon, lat = tmp$lat)
})

gaussianprofdf <- transform(gaussianprofdf,id=as.numeric(factor(juld)))

gaussdata <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- gaussiandf[gaussiandf$juld == gaussprofiles$juld[i],]
  tmp2 <- MLDdf[MLDdf$juld == gaussprofiles$juld[i],]
  data.frame(Fsurf = tmp$Fsurf, Zdemi = tmp$Zdemi, Fmax = tmp$Fmax,
             Zmax = tmp$Zmax, dz = tmp$dz, juld = tmp$juld, file = tmp$file,
             lon = tmp$lon, lat = tmp$lat, day = tmp$day, month = tmp$month,
             year = tmp$year, DOY = tmp$DOY, MLD = tmp2$MLD)
})

gaussdata <- transform(gaussdata,id=as.numeric(factor(juld)))

#SHAPE TEST FROM LAVIGNE 2015
HSCdf <- ldply(as.list(1:length(unique(gaussianprofdf$juld))), function(i){
  
  tmp <- gaussianprofdf[gaussianprofdf$id == i,]
  
  check_visu <- c(7,20,25,60,69,70,75,82,131,180,197,201,203,204,343,359,360,
                  361,364,369,371,373,379,390,514,515)
  check_visu <- 0
  minval <- min(tmp$fluo, na.rm = T)
  mindepth <- round(tmp$depth[which(tmp$fluo == minval)[1]])
  a <- seq(from = 0, to = mindepth, by = 10)#to take into account the increase in depth
  vec <- vector()
  NonNAindex <- min(which(!is.na(tmp$fluo)))
  
  for (j in 1:a[length(a)]/10){
    vec[j] <- mean(tmp$fluo[which(tmp$depth >= a[j] & tmp$depth <= a[j+1])], na.rm = T)
  }
  if(all(diff(vec) <= 0, na.rm = T)){
    result <- "HSC"
  }else if(max(tmp$fluo, na.rm = T) >= 2*tmp$fluo[which(!is.na(tmp$fluo))[1]]){
    result <- "DCM"
  }else if(max(tmp$fluo, na.rm = T) == tmp$fluo[NonNAindex]){
    result <- "HSC"
  }else if(i %in% check_visu ){
    result <- "Other"
  }else if(!(i %in% check_visu)){
    result <- "Modified_DCM"
  }
  
  data.frame(Category = result, juld = tmp$juld[1], file = tmp$Platform[1])
})

HSCdf <- transform(HSCdf,id=as.numeric(factor(juld)))


# GAUSSIAN VISU ONLY
#i <- 30
t1 <- which(HSCdf$Category == "HSC")
t2 <- which(HSCdf$Category == "DCM")
t3 <- which(HSCdf$Category == "Modified_DCM")
t4 <- which(HSCdf$Category == "Other")

#CHECK VISU
j <- 1
i <- t3[j]
tmp <- gaussianprofdf[gaussianprofdf$id==i,]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2, main=i, sub = 'test')
tmp2 <- approx(tab$x,tab$y, method="linear")
i <- i + 1
j <- j + 1

#TEMPORAL STATS ON PSEUDO GAUSSIAN PROFILES
temporal_stats <- cbind(gaussdata, HSCdf)
ggplot(temporal_stats, aes(month)) + geom_histogram(aes(fill=Category), 
                                                    binwidth = .5,
                                                    size = 1) +
  ylab("Counts") + scale_x_discrete("Month",
                                    limits = seq(1,12,by=1))


modified_DCM <- temporal_stats[temporal_stats$Category == "Modified_DCM",]
modified_DCM <- transform(modified_DCM,id=as.numeric(factor(juld)))
DCM <- temporal_stats[temporal_stats$Category == "DCM",]
DCM <- transform(DCM,id=as.numeric(factor(juld)))

#get DCM profiles
DCMprofiles <- ldply(as.list(1:length(DCM$juld)), function(i){
  tmp <- gaussprofiles[gaussprofiles$juld == DCM$juld[i],]
  tmp3 <- gaussdata[gaussdata$juld == DCM$juld[i],]
  tmp2 <- profiles[profiles$juld == DCM$juld[i],]
  n <- length(tmp2$depth)
  data.frame(depth = tmp2$depth, juld = tmp2$juld, fluo = tmp2$fluo, qc = tmp2$qc,
             day = tmp2$day, month = tmp2$month, year = tmp2$year, DOY = tmp2$DOY,
             Platform = tmp2$Platform, lon = tmp2$lon, lat = tmp2$lat, 
             Fsurf = rep(tmp3$Fsurf, n), Fmax = rep(tmp3$Fmax, n), Zmax = rep(tmp3$Zmax, n),
             dz = rep(tmp3$dz, n), Zdemi = rep(tmp3$Zdemi, n))
})

#get modified DCM
modified_DCMprofiles <- ldply(as.list(1:length(modified_DCM$juld)), function(i){
  tmp <- gaussprofiles[gaussprofiles$juld == modified_DCM$juld[i],]
  tmp3 <- gaussdata[gaussdata$juld == modified_DCM$juld[i],]
  tmp2 <- profiles[profiles$juld == modified_DCM$juld[i],]
  n <- length(tmp2$depth)
  data.frame(depth = tmp2$depth, juld = tmp2$juld, fluo = tmp2$fluo, qc = tmp2$qc,
             day = tmp2$day, month = tmp2$month, year = tmp2$year, DOY = tmp2$DOY,
             Platform = tmp2$Platform, lon = tmp2$lon, lat = tmp2$lat, 
             Fsurf = rep(tmp3$Fsurf, n), Fmax = rep(tmp3$Fmax, n), Zmax = rep(tmp3$Zmax, n),
             dz = rep(tmp3$dz, n), Zdemi = rep(tmp3$Zdemi, n))
})

DCMprofiles <- transform(DCMprofiles,id=as.numeric(factor(juld)))
modified_DCMprofiles <- transform(modified_DCMprofiles, id=as.numeric(factor(juld)))

# #Check Visu DCM
# tmp <- DCMprofiles[DCMprofiles$id==i,]#gappy data, no interpolation
# depthindex <- which.min(tmp$depth <= 100)
# if(check_data == "FLUO"){
#   tmp$fluo <- normalize(tmp$fluo)
# }
# tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
# x <- tab$x
# plot(y~x, data=tab, type="l", lwd = 2)
# lines(x = tab$x, fgauss(x,tmp$Fsurf[i], tmp$Zdemi[i],
#                         tmp$Fmax[i], tmp$Zmax[i], tmp$dz[i]),
#       col=4, add=T, xlim=range(tab$x), lwd = 2)
# i <- i + 1
# 
# #Check Visu modified_DCM
# tmp <- modified_DCMprofiles[modified_DCMprofiles$id==i,]#gappy data, no interpolation
# depthindex <- which.min(tmp$depth <= 100)
# if(check_data == "FLUO"){
#   tmp$fluo <- normalize(tmp$fluo)
# }
# tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
# x <- tab$x
# plot(y~x, data=tab, type="l", lwd = 2)
# lines(x = tab$x, fgauss(x,tmp$Fsurf[i], tmp$Zdemi[i],
#                         tmp$Fmax[i], tmp$Zmax[i], tmp$dz[i]),
#       col=4, add=T, xlim=range(tab$x), lwd = 2)
# i <- i + 1

################################ SIGMA COMPUTATION FOR DCM AND MODIFIED DCM PROFILES NOW

#reorder and sort by filename for rho_DCM
DCM <- DCM[order(DCM$file),]
modified_DCM <- modified_DCM[order(modified_DCM$file),]

#datadf is either DCM or modified_DCM
sigma_dcm_from_depth_dcm <- function(file, borne_inf, borne_sup, datadf){
  
  file <- paste0(file,"_Mprof.nc")
  ncfile <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  cycle_number <- ncvar_get(ncfile,"CYCLE_NUMBER")
  cycle_number <- cycle_number[N_PROF]
  
  FloatInfo<-list(N_PROF=N_PROF,
                  N_LEVELS=N_LEVELS,
                  juld = juld,
                  pres=pres,
                  lon=lon,
                  lat=lat,
                  cycle_number = cycle_number)
 
  tempdf <- ExtractVar("TEMP",FloatInfo)
  tempdf$value[which(tempdf$qc == 4 | tempdf$qc == 3)] <- NA
  psaldf <- ExtractVar("PSAL",FloatInfo)
  psaldf$value[which(psaldf$qc == 4 | psaldf$qc == 3)] <- NA

  df <- ldply(as.list(borne_inf:borne_sup), function(i){
    
    depth_dcm <- datadf$Zmax[i]
    julian_day <- datadf$juld[i]
    
    psaltmp <- psaldf[psaldf$juld == julian_day,]
    temptmp<- tempdf[tempdf$juld == julian_day,]
    
    indexdepth <- which.min(abs(temptmp$depth - depth_dcm))
    psaltmp <- psaltmp[indexdepth,]
    temptmp <- temptmp[indexdepth,]
    
    if(check_sigma =="TEOS10"){
      #DENSITY ANOMALY BASED ON TEOS-10
      psal <- gsw_SA_from_SP(psaltmp$value,depth_dcm,psaltmp$lon,psaltmp$lat)
      temp <- gsw_CT_from_t(psal,temptmp$value,depth_dcm) 
      sigma <- gsw_sigma0(psal,temp)
    }else{
      sigma <- swSigmaTheta(psaltmp$value,temptmp$value,temptmp$depth)
    }

    data.frame(sigma_dcm = sigma, filename = datadf$file[i],
               lon = datadf$lon[i], lat = datadf$lat[i],
               day = datadf$day[i], month = datadf$month[i], year = datadf$year[i],
               DOY = datadf$DOY[i], juld = datadf$juld[i])
  })
}

#FIRST : DCM DATA
datadf <- DCM

index_dcm <- which(duplicated(datadf$file) == FALSE)

sigma_dcm1 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[1]], index_dcm[1], index_dcm[2]-1, datadf)
sigma_dcm2 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[2]], index_dcm[2], index_dcm[3]-1, datadf)
sigma_dcm3 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[3]], index_dcm[3], index_dcm[4]-1, datadf)
sigma_dcm4 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[4]], index_dcm[4], length(datadf$juld), datadf)
sigma_dcm <- rbind(sigma_dcm1, sigma_dcm2, sigma_dcm3, sigma_dcm4)

rhoNA <- which(is.na(sigma_dcm$sigma_dcm) == TRUE)
dcm_juld_to_remove <- sigma_dcm$juld[rhoNA]#for below
#remove NA's
sigma_dcm <- sigma_dcm[complete.cases(sigma_dcm),]
sigma_dcm <- transform(sigma_dcm,id=as.numeric(factor(juld)))
sigmaDCM <- sigma_dcm

#---------------------#

#SECOND: MODIFIED DCM DATA
datadf <- modified_DCM

index_dcm <- which(duplicated(datadf$file) == FALSE)

sigma_dcm1 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[1]], index_dcm[1], index_dcm[2]-1, datadf)
sigma_dcm2 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[2]], index_dcm[2], index_dcm[3]-1, datadf)
sigma_dcm3 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[3]], index_dcm[3], index_dcm[4]-1, datadf)
sigma_dcm4 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[4]], index_dcm[4], length(datadf$juld), datadf)
sigma_dcm <- rbind(sigma_dcm1, sigma_dcm2, sigma_dcm3, sigma_dcm4)

rhoNA <- which(is.na(sigma_dcm$sigma_dcm) == TRUE)
modified_dcm_juld_to_remove <- sigma_dcm$juld[rhoNA]#for below
#remove NA's
sigma_dcm <- sigma_dcm[complete.cases(sigma_dcm),]
sigma_dcm <- transform(sigma_dcm,id=as.numeric(factor(juld)))

sigma_modified_DCM <- sigma_dcm


### DENSITY ANALYSIS ### -----------------

var_average <- function(dfin){
  dfin <- transform(dfin,id=as.numeric(factor(juld)))
  
  df <- ldply(as.list(1:length(unique(dfin$id))), function(i){
    #i <- 1
    tmp <- dfin[dfin$id == i,]
    tmp <- tmp[order(tmp$depth),]
    tmp2 <- aggregate(value ~ depth, data = tmp, FUN = mean)
    #tmp <- match_df(tmp2, tmp, on="depth")
    duplicate <- which(duplicated(tmp$depth))
    if (length(duplicate) != 0){
      tmp <- tmp[-duplicate,]
    }
    #tmp <- tmp[-duplicate,]
    tmp$value <- tmp2$value
    #clean (depth <0 without a QC = 4)
    #we remove them because it can induce further issues when computing
    #the potential density anomaly
    tmp <- tmp[!tmp$depth < 0,]
    #i <- i + 1
    data.frame(tmp)
  })
}

#Fonction qui renvoie un data frame basé sur l'anomalie de densité verticale (y axis) et la FLUO
ExtractDensity <- function(chladf, FloatInfo){
  
  tempdf <- ExtractVar("TEMP",FloatInfo)
  tempdf <- var_average(tempdf)#for averaging depth duplicates for temperature
  tempdf$value[which(tempdf$qc == 4 | tempdf$qc == 3)] <- NA
  psaldf <- ExtractVar("PSAL",FloatInfo)
  psaldf <- var_average(psaldf)#for averaging depth duplicates for salinity
  psaldf$value[which(psaldf$qc == 4 | psaldf$qc == 3)] <- NA
  
  subtempdf2 <- subset(tempdf, select = c("value","depth","juld","lon","lat"))
  colnames(subtempdf2)[which(colnames(subtempdf2)=="value")]<-"TEMP"
  subpsaldf2 <- subset(psaldf, select = c("value","depth","juld"))
  colnames(subpsaldf2)[which(colnames(subpsaldf2)=="value")]<-"PSAL"
  subchladf2 <- subset(chladf, select = c("fluo","depth","juld","qc"))
  joindf <- join(subtempdf2,subpsaldf2, by = c("depth","juld"))
  subchladf2 <- match_df(subchladf2,joindf, on = c("depth", "juld"))
  finaldf <- join(subchladf2,joindf, by = c("depth","juld"))
  
  if(check_sigma =="TEOS10"){
    #DENSITY ANOMALY BASED ON TEOS-10
    psal <- gsw_SA_from_SP(finaldf$PSAL,finaldf$depth,finaldf$lon,finaldf$lat)
    temp <- gsw_CT_from_t(psal,finaldf$TEMP,finaldf$depth)
    sigma <- gsw_sigma0(psal,temp)
  }else{
    sigma <- swSigmaTheta(finaldf$PSAL,finaldf$TEMP,finaldf$depth)
  }

  sigmadf <- data.frame(depth       = finaldf$depth,
                        sigma         = sigma,  
                        fluo          = finaldf$fluo,
                        TEMP          = temp,
                        PSAL          = psal,
                        juld          = finaldf$juld,
                        day           = month.day.year(finaldf$juld,c(1,1,1950))$day,
                        month         = month.day.year(finaldf$juld,c(1,1,1950))$month,
                        year          = month.day.year(finaldf$juld,c(1,1,1950))$year,
                        lon           = finaldf$lon,
                        lat           = finaldf$lat)
}

sigma_DCM_modified_densityprofiledf<-ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
  #Dimensions
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  cycle_number <- ncvar_get(ncfile,"CYCLE_NUMBER")
  cycle_number <- cycle_number[N_PROF]
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  FloatInfo<-list(N_PROF=N_PROF,
                  N_LEVELS=N_LEVELS,
                  juld = juld,
                  pres=pres,
                  lon=lon,
                  lat=lat,
                  cycle_number = cycle_number
  )
  
  #chladf <- profiledf[profiledf$Platform == as.numeric(unique(id)),]
  # MODIFY HERE
  chladf <- modified_DCMprofiles[modified_DCMprofiles$Platform == as.numeric(unique(id)),]
  
  #Extraction de l'anomalie de densité potentielle
  densitydf <- ExtractDensity(chladf, FloatInfo)
  
  #Construction du data frame final a lieu ici
  data.frame(depth         = densitydf$depth,
             density       = densitydf$sigma,
             juld          = densitydf$juld,
             fluo          = densitydf$fluo,
             day           = month.day.year(densitydf$juld,c(1,1,1950))$day,
             month         = month.day.year(densitydf$juld,c(1,1,1950))$month,
             year          = month.day.year(densitydf$juld,c(1,1,1950))$year,
             DOY           = as.integer(strftime(as.Date(densitydf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon           = densitydf$lon,
             lat           = densitydf$lat,
             Platform      = as.numeric(unique(id)),
             type          = "Argo")
})

sigma_DCM_densityprofiledf<-ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
  #Dimensions
  N_PROF   <- ncol(ncvar_get(ncfile,"PRES"))
  N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))
  juld     <- ncvar_get(ncfile,"JULD")
  pres     <- ncvar_get(ncfile,"PRES")
  lon      <- ncvar_get(ncfile,"LONGITUDE")
  lat      <- ncvar_get(ncfile,"LATITUDE")
  cycle_number <- ncvar_get(ncfile,"CYCLE_NUMBER")
  cycle_number <- cycle_number[N_PROF]
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  FloatInfo<-list(N_PROF=N_PROF,
                  N_LEVELS=N_LEVELS,
                  juld = juld,
                  pres=pres,
                  lon=lon,
                  lat=lat,
                  cycle_number = cycle_number
  )
  

  chladf <- DCMprofiles[DCMprofiles$Platform == as.numeric(unique(id)),]
  
  #Extraction de l'anomalie de densité potentielle
  densitydf <- ExtractDensity(chladf, FloatInfo)
  
  #Construction du data frame final a lieu ici
  data.frame(depth         = densitydf$depth,
             density       = densitydf$sigma,
             juld          = densitydf$juld,
             fluo          = densitydf$fluo,
             day           = month.day.year(densitydf$juld,c(1,1,1950))$day,
             month         = month.day.year(densitydf$juld,c(1,1,1950))$month,
             year          = month.day.year(densitydf$juld,c(1,1,1950))$year,
             DOY           = as.integer(strftime(as.Date(densitydf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon           = densitydf$lon,
             lat           = densitydf$lat,
             Platform      = as.numeric(unique(id)),
             type          = "Argo")
})


#pour supprimer les doublons qui apparaissent
sigma_DCM_densityprofiledf <- unique(sigma_DCM_densityprofiledf)
sigma_DCM_modified_densityprofiledf <- unique(sigma_DCM_modified_densityprofiledf)
#Assign a profile number according to the (unique due to two decimals) julian day 
sigma_DCM_densityprofiledf <- transform(sigma_DCM_densityprofiledf,id=as.numeric(factor(juld)))
sigma_DCM_modified_densityprofiledf <- transform(sigma_DCM_modified_densityprofiledf,id=as.numeric(factor(juld)))

#REMOVE PROFILES FOR WHICH JULD HAS TO BE REMOVED BECAUSE RHO COULD NOT 
# BE COMPUTED WITH CERTAINTY
# CHECK THAT AFTER IT IS WELL 228 AND 183

#DCM cleaning
for (i in 1:length(dcm_juld_to_remove)){
  sigma_DCM_densityprofiledf <- sigma_DCM_densityprofiledf[!(sigma_DCM_densityprofiledf$juld == dcm_juld_to_remove[i]),]
  DCMprofiles <- DCMprofiles[!(DCMprofiles$juld == dcm_juld_to_remove[i]),]
  DCM <- DCM[!DCM$juld == dcm_juld_to_remove[i],]
}

sigma_DCM_densityprofiledf <- transform(sigma_DCM_densityprofiledf,id=as.numeric(factor(juld)))
DCMprofiles <- transform(DCMprofiles, id = as.numeric(factor(juld)))
DCM <- transform(DCM, id = as.numeric(factor(juld)))

#Modified DCM cleaning
for (i in 1:length(modified_dcm_juld_to_remove)){
  sigma_DCM_modified_densityprofiledf <- sigma_DCM_modified_densityprofiledf[!(sigma_DCM_modified_densityprofiledf$juld == modified_dcm_juld_to_remove[i]),]
  modified_DCMprofiles <- modified_DCMprofiles[!(modified_DCMprofiles$juld == modified_dcm_juld_to_remove[i]),]
  modified_DCM <- modified_DCM[!(modified_DCM$juld == modified_dcm_juld_to_remove[i]),]
}

sigma_DCM_modified_densityprofiledf <- transform(sigma_DCM_modified_densityprofiledf,id=as.numeric(factor(juld)))
modified_DCMprofiles <- transform(modified_DCMprofiles, id = as.numeric(factor(juld)))
modified_DCM <- transform(modified_DCM, id = as.numeric(factor(juld)))

TOTAL_DCM_PROFILES <- rbind(DCMprofiles, modified_DCMprofiles)
TOTAL_DCM_PROFILES <- transform(TOTAL_DCM_PROFILES,id=as.numeric(factor(juld)))
TOTAL_DCM <- rbind(DCM, modified_DCM)
TOTAL_DCM <- transform(TOTAL_DCM, id = as.numeric(factor(juld)))
TOTAL_DCM <- TOTAL_DCM[order(TOTAL_DCM$id),]
TOTAL_SIGMA_PROFILES <- rbind(sigma_DCM_densityprofiledf, sigma_DCM_modified_densityprofiledf)
TOTAL_SIGMA_PROFILES <- transform(TOTAL_SIGMA_PROFILES, id = as.numeric(factor(juld)))
TOTAL_SIGMA_PROFILES <- TOTAL_SIGMA_PROFILES[order(TOTAL_SIGMA_PROFILES$id),]
TOTAL_SIGMA <- rbind(sigmaDCM, sigma_modified_DCM)
TOTAL_SIGMA <- transform(TOTAL_SIGMA, id = as.numeric(factor(juld)))
TOTAL_SIGMA <- TOTAL_SIGMA[order(TOTAL_SIGMA$id),]


TOTAL_DCM_PROFILES <- ddply(TOTAL_DCM_PROFILES,~month, transform, season=1*(month %in% c(12,1,2))+
                              2*(month %in% c(3,4,5 ))+
                              3*(month %in% c(6,7,8 ))+
                              4*(month %in% c(9,10,11 )))

TOTAL_DCM <- ddply(TOTAL_DCM,~month, transform, season=1*(month %in% c(12,1,2))+
                     2*(month %in% c(3,4,5 ))+
                     3*(month %in% c(6,7,8 ))+
                     4*(month %in% c(9,10,11 )))

TOTAL_SIGMA_PROFILES <- ddply(TOTAL_SIGMA_PROFILES,~month, transform, season=1*(month %in% c(12,1,2))+
                              2*(month %in% c(3,4,5 ))+
                              3*(month %in% c(6,7,8 ))+
                              4*(month %in% c(9,10,11 )))

TOTAL_SIGMA <- ddply(TOTAL_SIGMA,~month, transform, season=1*(month %in% c(12,1,2))+
                     2*(month %in% c(3,4,5 ))+
                     3*(month %in% c(6,7,8 ))+
                     4*(month %in% c(9,10,11 )))

#TEXT FILES creation FOR DIVA

# DCM ONLY

subDCM <- subset(DCM, select = c("lon", "lat", "juld", "Zmax","DOY", "year",
                                 "month", "day", "file"))
sub_sigmaDCM <- subset(sigmaDCM, select = c("lon", "lat", "juld", "sigma_dcm","DOY", "year",
                                        "month", "day", "filename"))

dcm_sigma_depth <- cbind(subDCM, sub_sigmaDCM$rho_sigma)

# write.table(subDCM, file="TRUE_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_rhoDCM, file="TRUE_DCM_RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(dcm_rho_depth, file="TRUE_DCM_DEPTH&RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

# MODIFIED DCM ONLY

sub_modified_DCM <- subset(modified_DCM, select = c("lon", "lat", "juld", "Zmax","DOY", "year",
                                                    "month", "day", "file"))
sub_sigma_modified_DCM <- subset(sigma_modified_DCM, select = c("lon", "lat", "juld", "sigma_dcm","DOY", "year",
                                                            "month", "day", "filename"))

modified_dcm_sigma_depth <- cbind(sub_modified_DCM, sub_sigma_modified_DCM$sigma_dcm)

# write.table(sub_modified_DCM, file="MODIFIED_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_rho_modified_DCM, file="MODIFIED_DCM_RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(modified_dcm_rho_depth, file="MODIFIED_DCM_DEPTH&RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#  row.names = FALSE, col.names = FALSE)

# ALL

sub_TOTAL_DCM <- subset(TOTAL_DCM, select = c("lon", "lat", "juld", "Zmax","DOY", "year",
                                              "month", "day", "file"))

sub_TOTAL_SIGMA <- subset(TOTAL_SIGMA, select = c("lon", "lat", "juld", "sigma_dcm","DOY", "year",
                                              "month", "day", "filename"))

TOTAL_DCM_DEPTH_SIGMA <- cbind(sub_TOTAL_DCM, sub_TOTAL_SIGMA$sigma_dcm)

# write.table(sub_TOTAL_DCM, file="TOTAL_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_TOTAL_RHO, file="TOTAL_DCM_RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(TOTAL_DCM_DEPTH_RHO, file="TOTAL_DCM_DEPTH&RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

NORM_DCM <- normalize(DCMprofiles$fluo)
DCMprofiles2 <- cbind(DCMprofiles, NORM_DCM)

NORM_modified_DCM <- normalize(modified_DCMprofiles$fluo)
modified_DCMprofiles2 <- cbind(modified_DCMprofiles$fluo)

NORM_TOT_DCM <- normalize(TOTAL_DCM_PROFILES$fluo)
TOTAL_DCM_PROFILES2 <- cbind(TOTAL_DCM_PROFILES, NORM_TOT_DCM)

NORM_TOT_RHO <- normalize(TOTAL_RHO_PROFILES$FLUO_ADJUSTED)
TOTAL_RHO_PROFILES2 <- cbind(TOTAL_RHO_PROFILES, NORM_TOT_RHO)

#USE FOR DCM PROFILES ONLY (meaning no MODIFIED_DCM)
NORM_RHO_DCM <- normalize(rho_DCM_densityprofiledf$FLUO_ADJUSTED)
RHOprofiles2 <- cbind(rho_DCM_densityprofiledf, NORM_RHO_DCM)


