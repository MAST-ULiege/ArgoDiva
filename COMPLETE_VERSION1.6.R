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
library(ggmap)
library(ggalt)

#/!!!!!!!!!!!!!!!!!!!!!!!\     CHECK BEFORE USE     /!!!!!!!!!!!!!!!!!!!!\

#check_data <- "FLUO"
check_data <- "CHLA"
check_sigma <- "TEOS10"
check_mld <- "ARGO"

#FLUORESCENCE PARAMETERS AVAILABLE IN METADATA FILES FROM EACH FLOATS
# param_6900807 <- c(0.0072, 47) #SCALE_CHLA, DARK_CHLA
# param_6901866 <- c(0.0073, 49)
# param_7900591 <- c(0.0121, 48)
# param_7900592 <- c(0.0121, 50)
# fluo_param <- list(param_6900807, param_6901866, param_7900591,
#                    param_7900592)

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
    
    lvar_dir       <- ncvar_get(ncfile,"DIRECTION")
    lvar_direction <- llply(lvar_dir,function(dirstring){
      strsplit(dirstring,split="")
    })
    lvar_direction <- unlist(lvar_direction)
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
                 dir      = lvar_direction[iprof],
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
  
  chladf <- ExtractVar("CHLA", FloatInfo)
  if (file == "7900591_Mprof.nc" | file == "7900592_Mprof.nc"){
    cdomdf <- as.data.frame(rep(NA, length(chladf$value)))
    colnames(cdomdf) <- "value"
  }else{
    cdomdf <- ExtractVar("CDOM", FloatInfo)
  }

  # #fluorescence conversion, respectively to each float
  
  # CHLA = (FLUO_CHLA - DARK_CHLA) * SCALE_CHLA
  # => FLUO_CHLA = (CHLA / SCALE_CHLA) + DARK_CHLA
  
  # if (check_data == "FLUO"){
  #   if (file == "6900807_Mprof.nc"){
  #     chladf$value <- (chladf$value / fluo_param[[1]][1]) + fluo_param[[1]][2] 
  #   }
  #   else if (file == "6901866_Mprof.nc"){
  #     chladf$value <- (chladf$value / fluo_param[[2]][1]) + fluo_param[[2]][2] 
  #   }
  #   else if (file == "7900591_Mprof.nc"){
  #     chladf$value <- (chladf$value / fluo_param[[3]][1]) + fluo_param[[3]][2] 
  #   }
  #   else if (file == "7900592_Mprof.nc"){
  #     chladf$value <- (chladf$value / fluo_param[[4]][1]) + fluo_param[[4]][2] 
  #   }}
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  data.frame(depth    = chladf$depth,
             juld     = chladf$juld,
             fluo     = chladf$value,#fluo means FChla (like in Xing et al. 2017)
             cdom     = cdomdf$value,
             qc       = chladf$qc,
             day      = month.day.year(chladf$juld,c(1,1,1950))$day,
             month    = month.day.year(chladf$juld,c(1,1,1950))$month,
             year     = month.day.year(chladf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = chladf$lon,
             lat      = chladf$lat,
             dir      = chladf$dir,
             Platform = as.numeric(unique(id)),
             type     = "Argo")
})

#REMOVE BAD DATA (I.E. QC = 4 FOR NEGATIVE SPIKES, JUMPS, ETC.)
#DELETION OF DEPTH WITH BAD DATA
profiles <- profiles[-(which(profiles$qc == 4)),] 

#REMOVE DESCENT PROFILES
profiles <- profiles[-(which(profiles$dir == "D")),]
rownames(profiles) <- NULL

#CREATION OF PROFILE IDs & REORDER DATA FRAME ACCORDING TO IT
profiles <- transform(profiles,id=as.numeric(factor(juld)))
profiles <- profiles[order(profiles$id),]

#MEDIAN FILTER TO REMOVE POSITIVE SPIKES (THEY MAY RETAIN INFORMATION BUT STRONG SPIKES 
#WILL NOT BE RETAINED FOR THIS STUDY)
smoothed_fluo <- ldply(as.list(unique(profiles$id)), function(i){
  tmp <- profiles[profiles$id == i,]
  tmp <- mmed(tmp$fluo, 5)
  data.frame(smoothed_fluo = tmp)
})

smoothed_cdom <- ldply(as.list(unique(profiles$id)), function(i){
  tmp <- profiles[profiles$id == i,]
  if (tmp$Platform[1] == 6900807 | tmp$Platform[1] == 6901866){
    tmp <- mmed(tmp$cdom, 5)
  }else{
    tmp <- tmp$cdom
  }
  data.frame(smoothed_cdom = tmp)
})

#REPLACE FLUO WITH SMOOTHED FLUO
profiles[,3] <- smoothed_fluo
profiles[,4] <- smoothed_cdom

#Construction of an ARGO-only dataframe -> first guess value for Zmax [NLS fit]
argodf<- ddply(profiles,~juld,summarize,
                    qc = qc[which.max(fluo)],
                    depthmax = depth[which.max(fluo)],
                    maxvalue = fluo[which.max(fluo)],
                    depthmin = depth[which.max(fluo):length(fluo)][which.min(fluo[which.max(fluo):length(fluo)])],
                    integration = sum(fluo),
                    bottomdepth = max(depth),
                    min_depth = min(depth),
                    cdom_sensor = is.na(cdom[1]),
                    dir = dir[1],
                    lon=mean(lon),
                    lat=mean(lat),
                    Platform = Platform[1],
                    day = day[1], month = month[1],
                    year = year[1], DOY = DOY[1])

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
             dir   = tempdf$dir,
             lon   = tempdf$lon,
             lat   = tempdf$lon,
             day   = month.day.year(tempdf$juld,c(1,1,1950))$day,
             month = month.day.year(tempdf$juld,c(1,1,1950))$month,
             year  = month.day.year(tempdf$juld,c(1,1,1950))$year,
             DOY   = as.integer(strftime(as.Date(tempdf$juld,origin = '1950-01-01'), 
                                         format ="%j")))#Day Of Year
  
})

#REMOVE DESCENT PROFILES
density_profiles <- density_profiles[-(which(density_profiles$dir == "D")),]
rownames(density_profiles) <- NULL

#IDs
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))
density_profiles <- density_profiles[order(density_profiles$id),]

#Get profiles id that do not answer to some criteria
criteriadf <- argodf[(argodf$min_depth > 5),]#all stuck profiles except one

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
rownames(density_profiles) <- NULL
rownames(profiles) <- NULL
rownames(argodf) <- NULL
profiles <- transform(profiles,id=as.numeric(factor(juld)))
argodf <- transform(argodf,id=as.numeric(factor(juld)))
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))

#COMPUTE MLD
sigma_criteria <- 0.03
depth_ref <- 10

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
                                       format ="%j")))
})

#REMOVE BAD PROFILES (NO MLD DATA AVAILABLE)
bad_mld <- which(is.na(MLDdf$MLD) == TRUE)
for (i in 1:length(bad_mld)){
  profiles <- profiles[!(profiles$id == bad_mld[i]),]
  density_profiles <- density_profiles[!(density_profiles$id == bad_mld[i]),]
}
argodf <- argodf[-(which(is.na(MLDdf$MLD) == TRUE)),]
MLDdf <- MLDdf[-(which(is.na(MLDdf$MLD) == TRUE)),]

rownames(density_profiles) <- NULL
rownames(profiles) <- NULL
rownames(argodf) <- NULL
rownames(MLDdf) <- NULL
MLDdf <- transform(MLDdf,id=as.numeric(factor(juld)))
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))
argodf <- transform(argodf,id=as.numeric(factor(juld)))
profiles <- transform(profiles,id=as.numeric(factor(juld)))

#FDOM-BASED CORRECTION AND MINIMUM-OFFSET CORRECTION
FDOM_OR_MINIMUM_OFFSET_CORRECTION <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  tmp <- profiles[profiles$id == i,]
  MaxDepth <- argodf$bottomdepth[i]#Max depth of the profile
  TopDepth <- argodf$depthmin[i]#Apparent minimum 
  
  if(all(is.na(tmp$cdom)) == TRUE){#no cdom sensor hence minimum-offset correction procedure
    offset <- tmp$fluo[which(tmp$depth == TopDepth)]
    fluo_cor <- tmp$fluo - offset
    fluo_cor[which(tmp$depth == TopDepth):length(tmp$depth)] <- 0
    slope_fdom <- NA
    C <- NA
  }else{#FDOM-based method
    calibrange <- tmp[which(tmp$depth == TopDepth):length(tmp$depth),]
    linearMod <- lm(fluo ~ cdom, data = calibrange)
    slope_fdom <- coef(linearMod)[[2]]
    C <- coef(linearMod)[[1]]
    fluo_cor <- tmp$fluo - slope_fdom*tmp$cdom - C
  }
  data.frame(depth = tmp$depth, fluo_cor = fluo_cor, 
             slope_fdom = rep(slope_fdom, length(fluo_cor)),
                              C = rep(C, length(fluo_cor)), juld = tmp$juld[1])
})


#REPLACE CORRECTED VALUES
profiles[,3] <- FDOM_OR_MINIMUM_OFFSET_CORRECTION$fluo_cor

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

#QUENCHING CORRECTION
fluo_NPQ <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  tmp <- profiles[profiles$id == i,]
  correction <- quenching_correction(tmp$fluo, tmp$depth, MLDdf$MLD[i])
  data.frame(fluo_NPQ = correction)
})
profiles$fluo <- fluo_NPQ$fluo_NPQ

#FIT
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
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$Platform[1])
  
})

sigmoidf <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- profiles[profiles$id==i,]
  depthindex <- which.min(tmp$depth <= 100)
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
    Zdemi <- 30#added due to a bug
    s <- -0.1
    res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res, "error")){
    Zdemi <- 30
    s <- -0.01
    res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res, "error")){
    Zdemi <- 30
    s <- -1
    res <- tryCatch(nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tab, control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  data.frame(Fsurf = coef(res)["Fsurf"], Zdemi = coef(res)["Zdemi"], s=coef(res)["s"], id=i,
             juld=tmp$juld[1], file=tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$Platform[1])
  
})


# #CHECK VISU
# tmp <- profiles[profiles$id==i,]#gappy data, no interpolation
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
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
  maxindex <- which.max(tmp$fluo)
  
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
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$Platform[1])
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
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$platform[1])
})

#GRAPHICAL RECAP OF TYPES OF PROFILES
# i <- 73
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
# 
# depthindex <- which.min(tmp$depth <= 100)
# tmp <- tmp[1:depthindex,]
# # a <- ggplot(tmp, aes(x = fluo, y = depth)) + geom_path() + scale_y_reverse()
# test <- as.data.frame(fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
#        gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i]))
# colnames(test) <- "fit"
# test <- cbind(tmp$depth, test)
# d <- ggplot(tmp, aes(x = fluo, y = depth)) + geom_path() + scale_y_reverse()+
#   geom_path(data = test, aes(x = fit, y = tmp$depth), colour = "green") +
#   xlab("Fluorescence (RFU)") + ylab("Depth (m)")
# test2 <- as.data.frame(fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
#                                 sigmoidf$s[i]))
# colnames(test2) <- "fit"
# test2 <- cbind(tmp$depth, test2)
# e <- ggplot(tmp, aes(x = fluo, y = depth)) + geom_path() + scale_y_reverse()+
#   geom_path(data = test2, aes(x = fit, y = tmp$depth), colour = "red") +
#   xlab("Fluorescence (RFU)") + ylab("Depth (m)")
# #89,197,123,73,21 -> For types of profiles


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
  }else if(max(tmp$fluo, na.rm = T)*0.90 < tmp$fluo[NonNAindex]){
    result <- "Other"
  }
  else{
    result <- "Modified_DCM"
  }
  data.frame(Category = result, juld = tmp$juld[1], file = tmp$Platform[1])
})

HSCdf <- transform(HSCdf,id=as.numeric(factor(juld)))

# GAUSSIAN VISU ONLY
t1 <- which(HSCdf$Category == "HSC")
t2 <- which(HSCdf$Category == "DCM")
t3 <- which(HSCdf$Category == "Modified_DCM")
t4 <- which(HSCdf$Category == "Other")

#CHECK VISU
j <- 1
i <- t4[j]
tmp <- gaussianprofdf[gaussianprofdf$id==i,]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2, main=i, sub = tmp$month[1])
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


#MAP OF DATA
# total <- rbind(DCM, modified_DCM)
# 
# total <- ddply(total,~month, transform, Season=1*(month %in% c(12,1,2))+
#                        2*(month %in% c(3,4,5 ))+
#                        3*(month %in% c(6,7,8 ))+
#                        4*(month %in% c(9,10,11 )))
# #Black Sea coordinates
# bs <- c(26.5,40,43,46)
# myMap <-get_map(location=bs, source="google", crop=FALSE, maptype = "terrain")
# ggmap(myMap) +
#   geom_point(aes(x=lon, y=lat, color=factor(Season)),
#              data = total, alpha = 1) + #+facet_grid(file) + scale_color_manual(values=c("red", "green", "yellow","blue"))
# labs(x = "Longitude (°)", y = "Latitude (°)", color = "Season\n") +
#   scale_color_manual(labels = c("Winter", "Spring","Summer","Autumn"),
#                      values = c("orangered","#7CAE00","#00BFC4","#C77CFF")) +
#   #annotate("point", x = 29, y = 43.1, colour = "black") +
#   geom_point(aes(x = 29, y = 43.1), color="black", shape = 4, size = 3) + facet_grid(~Category)
# 
# 


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

#Check Visu DCM
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

#Check Visu modified_DCM
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

sigma_mod_dcm1 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[1]], index_dcm[1], index_dcm[2]-1, datadf)
sigma_mod_dcm2 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[2]], index_dcm[2], index_dcm[3]-1, datadf)
sigma_mod_dcm3 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[3]], index_dcm[3], index_dcm[4]-1, datadf)
sigma_mod_dcm4 <- sigma_dcm_from_depth_dcm(datadf$file[index_dcm[4]], index_dcm[4], length(datadf$juld), datadf)
sigma_mod_dcm <- rbind(sigma_mod_dcm1, sigma_mod_dcm2, sigma_mod_dcm3, sigma_mod_dcm4)

rhoNA <- which(is.na(sigma_mod_dcm$sigma_dcm) == TRUE)
modified_dcm_juld_to_remove <- sigma_mod_dcm$juld[rhoNA]#for below
#remove NA's
sigma_mod_dcm <- sigma_mod_dcm[complete.cases(sigma_mod_dcm),]
sigma_mod_dcm <- transform(sigma_mod_dcm,id=as.numeric(factor(juld)))

sigma_modified_DCM <- sigma_mod_dcm


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
  
  sigmadf <- data.frame(depth         = finaldf$depth,
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

dcm_sigma_depth <- cbind(subDCM, sub_sigmaDCM$sigma_dcm)

# write.table(subDCM, file="TRUE_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_sigmaDCM, file="TRUE_DCM_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(dcm_sigma_depth, file="TRUE_DCM_DEPTH&SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
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
# write.table(sub_sigma_modified_DCM, file="MODIFIED_DCM_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(modified_dcm_sigma_depth, file="MODIFIED_DCM_DEPTH&SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
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
# write.table(sub_TOTAL_SIGMA, file="TOTAL_DCM_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(TOTAL_DCM_DEPTH_SIGMA, file="TOTAL_DCM_DEPTH&SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

NORM_DCM <- normalize(DCMprofiles$fluo)
DCMprofiles2 <- cbind(DCMprofiles, NORM_DCM)

NORM_modified_DCM <- normalize(modified_DCMprofiles$fluo)
modified_DCMprofiles2 <- cbind(modified_DCMprofiles$fluo)

NORM_TOT_DCM <- normalize(TOTAL_DCM_PROFILES$fluo)
TOTAL_DCM_PROFILES2 <- cbind(TOTAL_DCM_PROFILES, NORM_TOT_DCM)

NORM_TOT_SIGMA <- normalize(TOTAL_SIGMA_PROFILES$fluo)
TOTAL_SIGMA_PROFILES2 <- cbind(TOTAL_SIGMA_PROFILES, NORM_TOT_SIGMA)

#USE FOR DCM PROFILES ONLY (meaning no MODIFIED_DCM)
NORM_SIGMA_DCM <- normalize(sigma_DCM_densityprofiledf$fluo)
SIGMAprofiles2 <- cbind(sigma_DCM_densityprofiledf, NORM_SIGMA_DCM)

tmp <- TOTAL_SIGMA_PROFILES2[TOTAL_SIGMA_PROFILES2$year == 2017,]
# tmp <- TOTAL_DCM_PROFILES2[TOTAL_DCM_PROFILES2$year == 2017,]
Per1 <- tmp[tmp$DOY <= 10,]
Per2 <- tmp[tmp$DOY >= 30 & tmp$DOY <= 40,]
Per3 <- tmp[tmp$DOY >= 60 & tmp$DOY <= 70,]
Per4 <- tmp[tmp$DOY >= 90 & tmp$DOY <= 100,]
Per5 <- tmp[tmp$DOY >= 120 & tmp$DOY <= 130,]
Per6 <- tmp[tmp$DOY >= 150 & tmp$DOY <= 160,]
Per7 <- tmp[tmp$DOY >= 180 & tmp$DOY <= 190,]
Per8 <- tmp[tmp$DOY >= 210 & tmp$DOY <= 220,]
Per9 <- tmp[tmp$DOY >= 240 & tmp$DOY <= 250,]
Per10 <- tmp[tmp$DOY >= 270 & tmp$DOY <= 280,]
Per11 <- tmp[tmp$DOY >= 300 & tmp$DOY <= 310,]
Per12 <- tmp[tmp$DOY >= 330 & tmp$DOY <= 366,]

Per1$group <- c("1/10Jan")
Per2$group <- c("30Jan/9Feb")
Per3$group <- c("1/11Mar")
Per4$group <- c("31Mar/10Apr")
Per5$group <- c("30Apr/10May")
Per6$group <- c("30May/9Jun")
Per7$group <- c("29Jun/9Jul")
Per8$group <- c("29Jul/8Aug")
Per9$group <- c("28Aug/7Sep")
Per10$group <- c("27Sep/7Oct")
Per11$group <- c("27Oct/6Nov")
Per12$group <- c("26Nov/31Dec")

perioddf <- rbind(Per1, Per2, Per3, Per4, Per5, Per6, Per7, Per8, Per9, Per10, Per11, Per12)

#ordering facet
perioddf$group = factor(perioddf$group, levels=c('1/10Jan','30Jan/9Feb','1/11Mar',
                                                 '31Mar/10Apr','30Apr/10May','30May/9Jun',
                                                 '29Jun/9Jul','29Jul/8Aug','28Aug/7Sep',
                                                 '27Sep/7Oct','27Oct/6Nov','26Nov/31Dec'))

mld_model <- read.table("/home/flo/R/TFE/tfe/reanalysis_V4/MLD_analysis.txt",
                        header=T, sep="")
mld_2017 <- mld_model[1,]
mld_2016 <- mld_model[2,]
mld_2015 <- mld_model[3,]
mld_2014 <- mld_model[4,]

ggplot(perioddf, aes(x=NORM_TOT_SIGMA, y=density, group=juld)) +
  geom_point(size = 0.5) + facet_grid(~group) +
  ylab("Potential density anomaly (kg/m³)") + xlab("Fluorescence (rfu)") +
  scale_y_reverse() + 
  geom_hline(yintercept=mld_2017$MLD_MIN_ARGO, linetype="dashed") +
  geom_hline(yintercept=mld_2017$MLD_MAX_ARGO, linetype="dashed") +
  annotate("rect", xmin=-Inf, xmax = Inf, ymin = mld_2017$MLD_MIN_ARGO,
           ymax = mld_2017$MLD_MAX_ARGO, fill = "yellow",
           alpha = .5, color = NA) + 
  geom_hline(yintercept = mld_2017$MLD_MAX_MODEL, 
             colour ="red")


########## PAR ISOLUMES COMPUTATIONS ##########

PAR_files <- c("6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

PAR_profiles <- ldply(as.list(PAR_files),function(file){
  
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
  
  pardf <- ExtractVar("DOWNWELLING_PAR", FloatInfo) 
  
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  data.frame(depth    = pardf$depth,
             juld     = pardf$juld,
             par     = pardf$value,
             qc       = pardf$qc,
             day      = month.day.year(pardf$juld,c(1,1,1950))$day,
             month    = month.day.year(pardf$juld,c(1,1,1950))$month,
             year     = month.day.year(pardf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(pardf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = pardf$lon,
             lat      = pardf$lat,
             Platform = as.numeric(unique(id)),
             type     = "Argo")
})

#YEAR ANALYSIS
YEARPAR <- 2017

year_par_profiles <- PAR_profiles[PAR_profiles$year == YEARPAR,]
year_par_profiles <- transform(year_par_profiles,id=as.numeric(factor(juld)))

isolumes <- ldply(as.list(1:length(unique(year_par_profiles$juld))), function(i){
  tmp <- year_par_profiles[year_par_profiles$id == i,]
  par_surf <- tmp$par[1]
  p10 <- tmp$depth[which.max(tmp$par <= par_surf/10)]
  p1 <- tmp$depth[which.max(tmp$par <= par_surf/100)]
  data.frame(par_surf = par_surf, iso10 = p10, iso1 = p1,
             DOY = tmp$DOY[1], juld = tmp$juld[1])
})


tmp <- TOTAL_DCM_PROFILES2[TOTAL_DCM_PROFILES2$year == YEARPAR,]

crit_depth <- 100

#depth only between 0 and 100m
tmp <- ddply(TOTAL_DCM_PROFILES2,~juld,summarize,
             depth = depth[which(depth <= crit_depth)],
             DOY = DOY[which(depth <= crit_depth)],
             fluo = fluo[which(depth <= crit_depth)],
             juld = juld[which(depth <= crit_depth)])

Per1 <- tmp[tmp$DOY <= 20,]
Per2 <- tmp[tmp$DOY >= 30 & tmp$DOY <= 40,]
Per3 <- tmp[tmp$DOY >= 60 & tmp$DOY <= 70,]
Per4 <- tmp[tmp$DOY >= 90 & tmp$DOY <= 100,]
Per5 <- tmp[tmp$DOY >= 120 & tmp$DOY <= 130,]
Per6 <- tmp[tmp$DOY >= 150 & tmp$DOY <= 160,]
Per7 <- tmp[tmp$DOY >= 180 & tmp$DOY <= 190,]
Per8 <- tmp[tmp$DOY >= 210 & tmp$DOY <= 220,]
Per9 <- tmp[tmp$DOY >= 240 & tmp$DOY <= 250,]
Per10 <- tmp[tmp$DOY >= 270 & tmp$DOY <= 280,]
Per11 <- tmp[tmp$DOY >= 300 & tmp$DOY <= 310,]
Per12 <- tmp[tmp$DOY >= 330 & tmp$DOY <= 366,]

Per1$group <- c("1/10Jan")
Per2$group <- c("30Jan/9Feb")
Per3$group <- c("1/11Mar")
Per4$group <- c("31Mar/10Apr")
Per5$group <- c("30Apr/10May")
Per6$group <- c("30May/9Jun")
Per7$group <- c("29Jun/9Jul")
Per8$group <- c("29Jul/8Aug")
Per9$group <- c("28Aug/7Sep")
Per10$group <- c("27Sep/7Oct")
Per11$group <- c("27Oct/6Nov")
Per12$group <- c("26Nov/31Dec")

perioddf <- rbind(Per1, Per2, Per3, Per4, Per5, Per6, Per7, Per8, Per9, Per10, Per11, Per12)

doy1 <- unique(Per1$DOY)
doy2 <- unique(Per2$DOY)
doy3 <- unique(Per3$DOY)
doy4 <- unique(Per4$DOY)
doy5 <- unique(Per5$DOY)
doy6 <- unique(Per6$DOY)
doy7 <- unique(Per7$DOY)
doy8 <- unique(Per8$DOY)
doy9 <- unique(Per9$DOY)
doy10 <- unique(Per10$DOY)
doy11 <- unique(Per11$DOY)
doy12 <- unique(Per12$DOY)

j1 <- unique(Per1$juld)
j2 <- unique(Per2$juld)
j3 <- unique(Per3$juld)
j4 <- unique(Per4$juld)
j5 <- unique(Per5$juld)
j6 <- unique(Per6$juld)
j7 <- unique(Per7$juld)
j8 <- unique(Per8$juld)
j9 <- unique(Per9$juld)
j10 <- unique(Per10$juld)
j11 <- unique(Per11$juld)
j12 <- unique(Per12$juld)

#depth of dcm for each j period
jmax1 <- TOTAL_DCM[TOTAL_DCM$juld %in% j1,]$Zmax
jmax2 <- TOTAL_DCM[TOTAL_DCM$juld %in% j2,]$Zmax
jmax3 <- TOTAL_DCM[TOTAL_DCM$juld %in% j3,]$Zmax
jmax4 <- TOTAL_DCM[TOTAL_DCM$juld %in% j4,]$Zmax
jmax5 <- TOTAL_DCM[TOTAL_DCM$juld %in% j5,]$Zmax
jmax6 <- TOTAL_DCM[TOTAL_DCM$juld %in% j6,]$Zmax
jmax7 <- TOTAL_DCM[TOTAL_DCM$juld %in% j7,]$Zmax
jmax8 <- TOTAL_DCM[TOTAL_DCM$juld %in% j8,]$Zmax
jmax9 <- TOTAL_DCM[TOTAL_DCM$juld %in% j9,]$Zmax
jmax10 <- TOTAL_DCM[TOTAL_DCM$juld %in% j10,]$Zmax
jmax11 <- TOTAL_DCM[TOTAL_DCM$juld %in% j11,]$Zmax
jmax12 <- TOTAL_DCM[TOTAL_DCM$juld %in% j12,]$Zmax

sj1 <- c(mean(jmax1), sd(jmax1))
sj2 <- c(mean(jmax2), sd(jmax2))
sj3 <- c(mean(jmax3), sd(jmax3))
sj4 <- c(mean(jmax4), sd(jmax4))
sj5 <- c(mean(jmax5), sd(jmax5))
sj6 <- c(mean(jmax6), sd(jmax6))
sj7 <- c(mean(jmax7), sd(jmax7))
sj8 <- c(mean(jmax8), sd(jmax8))
sj9 <- c(mean(jmax9), sd(jmax9))
sj10 <- c(mean(jmax10), sd(jmax10))
sj11 <- c(mean(jmax11), sd(jmax11))
sj12 <- c(mean(jmax12), sd(jmax12))

sjmean <- as.data.frame(c(sj1[1], sj2[1], sj3[1],sj4[1],sj5[1],sj6[1],
                          sj7[1],sj8[1],sj9[1],sj10[1],sj11[1],sj12[1]))
colnames(sjmean) <- "mean"

sjsd <- as.data.frame(c(sj1[2], sj2[2], sj3[2],sj4[2],sj5[2],sj6[2],
                        sj7[2],sj8[2],sj9[2],sj10[2],sj11[2],sj12[2]))

colnames(sjsd) <- "sd"

doy1PAR <- isolumes[isolumes$DOY %in% doy1,]
doy2PAR <- isolumes[isolumes$DOY %in% doy2,]
doy3PAR <- isolumes[isolumes$DOY %in% doy3,]
doy4PAR <- isolumes[isolumes$DOY %in% doy4,]
doy5PAR <- isolumes[isolumes$DOY %in% doy5,]
doy6PAR <- isolumes[isolumes$DOY %in% doy6,]
doy7PAR <- isolumes[isolumes$DOY %in% doy7,]
doy8PAR <- isolumes[isolumes$DOY %in% doy8,]
doy9PAR <- isolumes[isolumes$DOY %in% doy9,]
doy10PAR <- isolumes[isolumes$DOY %in% doy10,]
doy11PAR <- isolumes[isolumes$DOY %in% doy11,]
doy12PAR <- isolumes[isolumes$DOY %in% doy12,]

#compute mean iso 10% et 1% 
doy1meaniso <- c(mean(doy1PAR$iso10), mean(doy1PAR$iso1))
doy2meaniso <- c(mean(doy2PAR$iso10), mean(doy2PAR$iso1))
doy3meaniso <- c(mean(doy3PAR$iso10), mean(doy3PAR$iso1))
doy4meaniso <- c(mean(doy4PAR$iso10), mean(doy4PAR$iso1))
doy5meaniso <- c(mean(doy5PAR$iso10), mean(doy5PAR$iso1))
doy6meaniso <- c(mean(doy6PAR$iso10), mean(doy6PAR$iso1))
doy7meaniso <- c(mean(doy7PAR$iso10), mean(doy7PAR$iso1))
doy8meaniso <- c(mean(doy8PAR$iso10), mean(doy8PAR$iso1))
doy9meaniso <- c(mean(doy9PAR$iso10), mean(doy9PAR$iso1))
doy10meaniso <- c(mean(doy10PAR$iso10), mean(doy10PAR$iso1))
doy11meaniso <- c(mean(doy11PAR$iso10), mean(doy11PAR$iso1))
doy12meaniso <- c(mean(doy12PAR$iso10), mean(doy12PAR$iso1))

NORM_TOT_DCM_DEPTH <- normalize(perioddf$fluo)
perioddf2 <- cbind(perioddf, NORM_TOT_DCM_DEPTH)

index_group <- which(duplicated(perioddf$group) == FALSE)

alliso10 <- as.data.frame(c(doy1meaniso[1], doy2meaniso[1],doy3meaniso[1],
                            doy4meaniso[1],doy5meaniso[1],doy6meaniso[1],
                            doy7meaniso[1],doy8meaniso[1],doy9meaniso[1],
                            doy10meaniso[1],doy11meaniso[1],doy12meaniso[1]))
colnames(alliso10) <- "iso10"
alliso1 <- as.data.frame(c(doy1meaniso[2], doy2meaniso[2],doy3meaniso[2],
                           doy4meaniso[2],doy5meaniso[2],doy6meaniso[2],
                           doy7meaniso[2],doy8meaniso[2],doy9meaniso[2],
                           doy10meaniso[2],doy11meaniso[2],doy12meaniso[2]))
colnames(alliso1) <- "iso1"

#ordering facet
perioddf2$group = factor(perioddf2$group, levels=c('1/10Jan','30Jan/9Feb','1/11Mar',
                                                   '31Mar/10Apr','30Apr/10May','30May/9Jun',
                                                   '29Jun/9Jul','29Jul/8Aug','28Aug/7Sep',
                                                   '27Sep/7Oct','27Oct/6Nov','26Nov/31Dec'))

ggplot(perioddf2, aes(x=NORM_TOT_DCM_DEPTH, y=depth, group=juld)) +
  geom_path(size = 0.5) + facet_grid(~group) +
  ylab("Depth (m)") + xlab("Fluorescence (rfu)") +
  scale_y_reverse() + 
  geom_hline(data = data.frame(yint=as.vector(alliso10$iso10), 
                               group=c('1/10Jan','30Jan/9Feb','1/11Mar',
                                       '31Mar/10Apr','30Apr/10May','30May/9Jun',
                                       '29Jun/9Jul','29Jul/8Aug','28Aug/7Sep',
                                       '27Sep/7Oct','27Oct/6Nov','26Nov/31Dec')),
             aes(yintercept=yint), colour = "red") +
  geom_hline(data = data.frame(yint=as.vector(alliso1$iso1), 
                               group=c('1/10Jan','30Jan/9Feb','1/11Mar',
                                       '31Mar/10Apr','30Apr/10May','30May/9Jun',
                                       '29Jun/9Jul','29Jul/8Aug','28Aug/7Sep',
                                       '27Sep/7Oct','27Oct/6Nov','26Nov/31Dec')),
             aes(yintercept=yint), colour = "red", linetype = "dotted") +
  geom_hline(data = data.frame(yint=as.vector(sjmean$mean), 
                               group=c('1/10Jan','30Jan/9Feb','1/11Mar',
                                       '31Mar/10Apr','30Apr/10May','30May/9Jun',
                                       '29Jun/9Jul','29Jul/8Aug','28Aug/7Sep',
                                       '27Sep/7Oct','27Oct/6Nov','26Nov/31Dec')),
             aes(yintercept=yint), colour = "yellow") 


#Add l'incertitude sur la moyenne des profondeurs
#Régler le problème de fin Novembre
#Refaire les figures 3b de Navarro
#Clôturer les fichiers !! -> Lancer les analyses DIVA



