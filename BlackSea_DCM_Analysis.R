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
library(RColorBrewer)
#library(pbapply) #can be useful for time evolution of processing functions

#FLUORESCENCE PARAMETERS AVAILABLE IN METADATA FILES FROM EACH FLOATS (DARK_CHLA and SCALE_CHLA)

#For the definition of potential density anomaly
check_sigma <- "TEOS10"

#CHLA-equipped ARGO floats
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

#Extract variable
ExtractVar <- function(Var,FloatInfo){
  with(FloatInfo,{
    # This function should return a dataframe for variable with value, qc, iprofile and ilevel
    lvar             <- ncvar_get(ncfile,Var)
    lvar_qc          <- ncvar_get(ncfile,paste0(Var,"_QC"))
    lvar_qctab       <- llply(lvar_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab <- do.call(cbind,lvar_qctab)
    
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

#5-points median filter
mmed <- function(x,n=5){runmed(x,n)}

#NORMALIZATION (FOR FLUO DATA BEFORE GAUSSIAN FITTING FOR INSTANCE)
normalize <- function(data){
  data <- (data - min(data, na.rm = T))/(max(data, na.rm = T) - min(data, na.rm = T))
}

#Get each profile from ARGO data 
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
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  data.frame(depth    = chladf$depth,
             juld     = chladf$juld,
             fluo     = chladf$value,#fluo means FChla (Xing et al. 2017)
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
rownames(profiles) <- NULL#lines index can be problematic when deleting indexes and searching for max, min, etc. AFTER

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

#EXAMPLE DCM
# tmp <- profiles[profiles$id == 110,]
# ggplot(tmp, aes(x = fluo, y = depth)) + geom_path() + scale_y_reverse() +
#   xlab(expression(Chlorophyll~a~(mg/m^3))) + ylab("Depth (m)") + 
#   theme(text=element_text(size=12)) + geom_vline(xintercept = 0, colour = "red",
#                                                  lty = "dashed")

#Construction of an ARGO-only dataframe -> first guess value for Zmax [NLS fit]
argodf<- ddply(profiles,~juld,summarize,
               qc = qc[which.max(fluo)],
               depthmax = depth[which.max(fluo)],
               maxvalue = fluo[which.max(fluo)],
               depthmin = depth[which.max(fluo):length(fluo)][which.min(fluo[which.max(fluo):length(fluo)])],#to not take into account surface depth where chlorophyll can be minimum
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
    sigma <- swSigmaTheta(psaldf$value,tempdf$value,tempdf$depth)#Other definition of sigma_theta
  }
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER") 
  
  data.frame(sigma = sigma,
             depth = tempdf$depth, 
             juld  = tempdf$juld, 
             dir   = tempdf$dir,
             lon   = tempdf$lon,
             lat   = tempdf$lat,
             day   = month.day.year(tempdf$juld,c(1,1,1950))$day,
             month = month.day.year(tempdf$juld,c(1,1,1950))$month,
             year  = month.day.year(tempdf$juld,c(1,1,1950))$year,
             DOY   = as.integer(strftime(as.Date(tempdf$juld,origin = '1950-01-01'), 
                                         format ="%j")),
             Platform = as.numeric(unique(id)))#Day Of Year
  
})

#REMOVE DESCENT PROFILES
density_profiles <- density_profiles[-(which(density_profiles$dir == "D")),]
rownames(density_profiles) <- NULL

#IDs
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))
density_profiles <- density_profiles[order(density_profiles$id),]

#Get profiles id that do not answer to some criteria
criteriadf <- argodf[(argodf$min_depth > 5),]#all stuck profiles except one (add QC on pressure to get rid of them directly)

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

#COMPUTE MLD (criteria can thus be modified according to the definition of MLD that you choose)
sigma_criteria <- 0.03
depth_ref <- 10

MLDdf <- ldply(as.list(1:length(unique(density_profiles$id))), function(i){
  
  tmp <- density_profiles[density_profiles$id == i,]
  rownames(tmp) <- NULL
  
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
  
  if(is.na(MLD) == TRUE){
    sigmaMaxMLD <- NA
  }else{
    sigmaMaxMLD <- mean(tmp$sigma[which(tmp$depth == MLD)], na.rm =T)
  }
  
  show(i)
  
  data.frame(sigma_surface = sigma_surface, MLD = MLD, sigmaMaxMLD = sigmaMaxMLD,
             juld = tmp$juld[1], day = month.day.year(tmp$juld[1],c(1,1,1950))$day,
             month = month.day.year(tmp$juld[1],c(1,1,1950))$month,
             year = month.day.year(tmp$juld[1],c(1,1,1950))$year,
             DOY = as.integer(strftime(as.Date(tmp$juld[1],origin = '1950-01-01'), 
                                       format ="%j")),
             Platform = tmp$Platform[1])
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

#QUENCHING CORRECTION
#NOTE : This function was given to me by Marin Cornec (PhD Student at LOV)
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

#FDOM CORRECTION EXAMPLE
# tmp <- profiles[profiles$id == 150,]
# a <- ggplot(tmp, aes(x = fluo, y  = depth)) + geom_point(colour = "red") + 
#   scale_y_reverse() + geom_vline(xintercept = 0, colour = "black", lty = "dashed") +
#   geom_point(data = tmp2, aes(x = fluo, y = depth)) +
#   xlab(expression(Chlorophyll~a~(mg/m^3))) + ylab("Depth (m)") + 
#   theme(text=element_text(size=12)) + 
#   geom_hline(yintercept = 70.7, colour = "black", lty = "dashed")
#   
# b <- ggplot(tmp, aes(x = cdom, y  = depth)) + geom_point(colour = "purple") + 
#   scale_y_reverse() + geom_hline(yintercept = 70.7, colour = "black", lty = "dashed") +
#   xlab(expression(CDOM~(ppb))) + ylab("Depth (m)") + 
#   theme(text=element_text(size=12)) +
#   theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
#         axis.ticks.y=element_blank()) 
# 
# grid.arrange(a,b,nrow = 1, ncol=2)

#QUENCHING CORRECTION APPLICATION
fluo_NPQ <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  tmp <- profiles[profiles$id == i,]
  correction <- quenching_correction(tmp$fluo, tmp$depth, MLDdf$MLD[i])
  data.frame(fluo_NPQ = correction)
})
profiles$fluo <- fluo_NPQ$fluo_NPQ

#NPQ CORRECTION EXAMPLE
# tmp <- profiles[profiles$id ==6,]
# a <- ggplot(tmp2, aes(x = fluo, y = depth)) + geom_path() + 
#   scale_y_reverse() + ylim(c(100,0)) + 
#   xlab(expression(Chlorophyll~a~(mg/m^3))) + ylab("Depth (m)") +
#   theme(text=element_text(size=12)) 
# 
# b <- ggplot(tmp, aes(x = fluo, y = depth)) + geom_path() + 
#   scale_y_reverse() + ylim(c(100,0)) + 
#   xlab(expression(Chlorophyll~a~(mg/m^3))) + ylab("Depth (m)") +
#   theme(text=element_text(size=12)) +
#   theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# grid.arrange(a,b,nrow = 1, ncol = 2)

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

#COEFFICIENT OF DETERMINATION AFTER FIT
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
R_crit <- 0.8 #80% OF EXPLAINED VARIANCE (90% FOR NAVARRO ET AL. 2013)

classification <- ldply(as.list(1:length(unique(profiles$id))), function(i){
  
  tmp <- gaussiandf[gaussiandf$id==i,]
  
  if(Rcoefdf$rcoef_sigmoid[i] & Rcoefdf$rcoef_gaussian[i] < R_crit){
    classif <- "other"
    score <- "NA"
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

#Count profiles
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

#reorder by juld
gaussianprofdf <- transform(gaussianprofdf,id=as.numeric(factor(juld)))

#get gaussian parameters
gaussdata <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- gaussiandf[gaussiandf$juld == gaussprofiles$juld[i],]
  tmp2 <- MLDdf[MLDdf$juld == gaussprofiles$juld[i],]
  data.frame(Fsurf = tmp$Fsurf, Zdemi = tmp$Zdemi, Fmax = tmp$Fmax,
             Zmax = tmp$Zmax, dz = tmp$dz, juld = tmp$juld, file = tmp$file,
             lon = tmp$lon, lat = tmp$lat, day = tmp$day, month = tmp$month,
             year = tmp$year, DOY = tmp$DOY, MLD = tmp2$MLD)
})

#reorder
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
    result <- "Undetermined"
  }
  else{
    result <- "Modified DCM"
  }
  data.frame(Category = result, juld = tmp$juld[1], file = tmp$Platform[1])
})

HSCdf <- transform(HSCdf,id=as.numeric(factor(juld)))

# GAUSSIAN VISU ONLY
t1 <- which(HSCdf$Category == "HSC")
t2 <- which(HSCdf$Category == "DCM")
t3 <- which(HSCdf$Category == "Modified DCM")
t4 <- which(HSCdf$Category == "Undetermined")

#TEMPORAL STATS ON PSEUDO GAUSSIAN PROFILES
temporal_stats <- cbind(gaussdata, HSCdf)
# 
# ggplot(temporal_stats, aes(month)) + geom_histogram(aes(fill=Category), 
#                                                     binwidth = .5,
#                                                     size = 1,
#                                                     col = "black") +
#   ylab("Number of profiles") + scale_x_discrete("Month",
#                                                 limits = seq(1,12,by=1)) +
#   theme(text = element_text(size=16)) + 
#   scale_fill_manual(values = c("red", "red4","steelblue2","steelblue4"), labels = c("HSC  ",
#                                                                                     "Undetermined   ",
#                                                                                     "Modified DCM   ",
#                                                                                     "DCM")) +
#   theme(legend.position="bottom") + theme(legend.title=element_blank()) 
# 


#TEMPORAL STATS ON ALL PROFILES
# ggplot(classification, aes(month)) + geom_histogram(aes(fill=classif), 
#                                                     binwidth = .5,
#                                                     size = 1,
#                                                     col = "black") +
#   ylab("Number of profiles") + scale_x_discrete("Month",
#                                                 limits = seq(1,12,by=1)) +
#   theme(text = element_text(size=16)) + 
#   scale_fill_manual(values = c("red","steelblue2","steelblue4"),
#                     labels = c("Sigmoid  ", "Gaussian  ", "Other")) +
#   theme(legend.position="bottom") + theme(legend.title=element_blank())


modified_DCM <- temporal_stats[temporal_stats$Category == "Modified DCM",]
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

################################ SIGMA COMPUTATION FOR DCM AND MODIFIED DCM PROFILES NOW

#reorder and sort by filename for rho_DCM
DCM <- DCM[order(DCM$file),]
modified_DCM <- modified_DCM[order(modified_DCM$file),]

#get sigma dcm from depth dcm
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

#////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\\\\\\\\\\\\\\\\\\\
#BAD CODE BUT WORKING HERE
#SPLITTING DUE TO A BUG WITH CHANGE OF FILENAME THUS THIS IS NOT APPLICABLE IF WE HAVE MORE THAN 4 ARGO FILES
#///////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\\\\\\\\\\\\\\\\\\\
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
    tmp <- dfin[dfin$id == i,]
    tmp <- tmp[order(tmp$depth),]
    tmp2 <- aggregate(value ~ depth, data = tmp, FUN = mean)
    duplicate <- which(duplicated(tmp$depth))
    if (length(duplicate) != 0){
      tmp <- tmp[-duplicate,]
    }
    tmp$value <- tmp2$value
    #clean (depth <0 without a QC = 4)
    #we remove them because it can induce further issues when computing
    #the potential density anomaly
    tmp <- tmp[!tmp$depth < 0,]
    #i <- i + 1
    data.frame(tmp)
  })
}

#Function producing a data frame based on potential density anomaly and fluo
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
  
  densitydf <- ExtractDensity(chladf, FloatInfo)
  
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


#Remove exact duplicates
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

sub_sigmaDCM$Category <- "DCM"

dcm_sigma_depth <- cbind(subDCM, sub_sigmaDCM$sigma_dcm)

# write.table(subDCM, file="TRUE_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_sigmaDCM, file="TRUE_DCM_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(dcm_sigma_depth, file="TRUE_DCM_DEPTH_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

# MODIFIED DCM ONLY

sub_modified_DCM <- subset(modified_DCM, select = c("lon", "lat", "juld", "Zmax","DOY", "year",
                                                    "month", "day", "file"))
sub_sigma_modified_DCM <- subset(sigma_modified_DCM, select = c("lon", "lat", "juld", "sigma_dcm","DOY", "year",
                                                                "month", "day", "filename"))

sub_sigma_modified_DCM$Category <- "MODDCM"

modified_dcm_sigma_depth <- cbind(sub_modified_DCM, sub_sigma_modified_DCM$sigma_dcm)

# write.table(sub_modified_DCM, file="MODIFIED_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_sigma_modified_DCM, file="MODIFIED_DCM_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(modified_dcm_sigma_depth, file="MODIFIED_DCM_DEPTH_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#  row.names = FALSE, col.names = FALSE)

# ALL

sub_TOTAL_DCM <- subset(TOTAL_DCM, select = c("lon", "lat", "juld", "Zmax","DOY", "year",
                                              "month", "day", "file"))

sub_TOTAL_SIGMA <- subset(TOTAL_SIGMA, select = c("lon", "lat", "juld", "sigma_dcm","DOY", "year",
                                                  "month", "day", "filename"))

TOTAL_DCM_DEPTH_SIGMA <- cbind(sub_TOTAL_DCM, sub_TOTAL_SIGMA$sigma_dcm)

#IN FACT, TOTAL2 NOT NEEDED (SEE TOTAL_DCM_DEPTH_SIGMA) UNLESS YOU WANT THE CATEGORY SO KEEP IT 
TOTAL2 <- rbind(sub_sigmaDCM, sub_sigma_modified_DCM)
TOTAL2 <- TOTAL2[order(TOTAL2$juld),]
TOTAL_DCM_DEPTH_SIGMA <- TOTAL_DCM_DEPTH_SIGMA[order(TOTAL_DCM_DEPTH_SIGMA$juld),]
TOTAL2$Zmax <- TOTAL_DCM_DEPTH_SIGMA$Zmax
TOTAL2 <- transform(TOTAL2, id = as.numeric(factor(juld)))

#ADD SIGMA MLD PER PROFILE ("pseudo" MAX MLD)
MLDdf2 <- MLDdf[MLDdf$juld %in% TOTAL2$juld,]
MLDdf2 <- transform(MLDdf2, id = as.numeric(factor(juld)))
TOTAL2$sigma_MLD <- MLDdf2$sigmaMaxMLD
range_jan_feb2014 <- range(MLDdf[MLDdf$year == 2014 & MLDdf$month %in% c(1,2),]$sigmaMaxMLD)
range_jan_feb2017 <- range(MLDdf[MLDdf$year == 2017 & MLDdf$month %in% c(1,2),]$sigmaMaxMLD)
range_jan_feb2015 <- range(MLDdf[MLDdf$year == 2015 & MLDdf$month %in% c(1,2),]$sigmaMaxMLD)
range_jan_feb2016 <- range(MLDdf[MLDdf$year == 2016 & MLDdf$month %in% c(1,2),]$sigmaMaxMLD)
rangedf <- data.frame(year = c(2014,2015,2016,2017), inf = c(range_jan_feb2014[1],
                                                             range_jan_feb2015[1],range_jan_feb2016[1],range_jan_feb2017[1]), sup =
                        c(range_jan_feb2014[2],
                          range_jan_feb2015[2],range_jan_feb2016[2],range_jan_feb2017[2]))

#MLD each platform and each year (DCM OR NOT HENCE 622 PROFILES SPANNING 2014 TO 2017)
MLDdf3 <- ddply(MLDdf, ~year~Platform, summarize,
                max = max(MLD),
                sigmaMaxMLD = sigmaMaxMLD[which.max(MLD)])

#Ratio sigmaDCM/sigamMAXMLD (per year and per platform)
sigma_ratio <- ldply(as.list(1:length(MLDdf2$id)), function(i){
  tmp <- MLDdf2[MLDdf2$id == i,]
  tmp2 <- TOTAL2[TOTAL2$id == i,]
  sigmaMaxMLD <- MLDdf3$sigmaMaxMLD[which(MLDdf3$year == tmp$year & MLDdf3$Platform == tmp$Platform)]
  ratio = tmp2$sigma_dcm/sigmaMaxMLD
  data.frame(ratio = ratio, juld = tmp$juld[1])
})

TOTAL2$sigma_ratio <- sigma_ratio$ratio

ggplot(TOTAL2, aes(sigma_ratio, fill=Category)) +  geom_histogram(data = TOTAL2, col = "black")+ 
  theme(legend.position="none") +
  ylab("Number of profiles") + xlab(expression(paste(sigma["DCM"],"/",sigma["MLD-MAX"]))) + 
  scale_fill_manual(values=c("steelblue4","steelblue2")) +
  theme(text = element_text(size=13)) 

#REARRANGE TOTAL2 AND PUT DEPTH IN IT THEN CHECK WITH TOTAL_DCM_DETH_SIGMA

# write.table(sub_TOTAL_DCM, file="TOTAL_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_TOTAL_SIGMA, file="TOTAL_DCM_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(TOTAL_DCM_DEPTH_SIGMA, file="TOTAL_DCM_DEPTH_SIGMA.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)


TOTAL2 <- ddply(TOTAL2,~month, transform, Season=1*(month %in% c(12,1,2)) +
                  2*(month %in% c(3,4,5 )) +
                  3*(month %in% c(6,7,8 )) +
                  4*(month %in% c(9,10,11 )))

# write.table(TOTAL2, file="FINAL_CHECK.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 

#MLD time of event
# new_time <- as.Date(MLDdf2$juld, origin = as.Date("1950-01-01"))
# TOTAL2$date <- new_time
# TOTAL2$MLD <- MLDdf2$MLD
# 
# navarro8 <- ddply(TOTAL2, ~month~year, summarize,
#                   mld = mean(MLD),
#                   dcm = mean(Zmax))
# 
# ggplot(TOTAL2, aes(x = date, y = MLD)) +
#   geom_path() + geom_path(data = TOTAL2, aes(x = date, y = Zmax), colour = "red")


# chloro_dcm <- ddply(TOTAL_DCM_PROFILES, ~month, summarize,
#                     max = max(fluo),
#                     month = unique(month))

# write.table(TOTAL2, file="FULL_TOTAL_RATIO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
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

tmp <- TOTAL_SIGMA_PROFILES2[TOTAL_SIGMA_PROFILES2$year == 2017,]#MODIFY YEAR HERE

Per2 <- tmp[tmp$month == 2,]
Per3 <- tmp[tmp$month == 3,]
Per4 <- tmp[tmp$month == 4,]
Per5 <- tmp[tmp$month == 5,]
Per6 <- tmp[tmp$month == 6,]
Per7 <- tmp[tmp$month == 7,]
Per8 <- tmp[tmp$month == 8,]
Per9 <- tmp[tmp$month == 9,]

Per2$group <- c("February")
Per3$group <- c("March")
Per4$group <- c("April")
Per5$group <- c("May")
Per6$group <- c("June")
Per7$group <- c("July")
Per8$group <- c("August")
Per9$group <- c("September")

perioddf <- rbind(Per2, Per3, Per4, Per5, Per6, Per7, Per8, Per9)

#ordering facet
perioddf$group = factor(perioddf$group, 
                        levels=c("February", "March", "April","May","June","July","August","September"))

NORM_TOT_SIGMABIS <- normalize(perioddf$fluo)

perioddf2 <- cbind(perioddf, NORM_TOT_SIGMABIS)

info2017 <- ldply(as.list(2:9), function(i){
  tmp2 <- tmp[tmp$month == i,]
  tmp2 <- transform(tmp2,id=as.numeric(factor(juld)))
  sigma_max <- ddply(tmp2, ~juld, summarize,
                     maxfluo = mean(density[which.max(fluo)], na.rm = T))
  mean <- mean(sigma_max$maxfluo)
  sd <- sd(sigma_max$maxfluo)
  data.frame(mean = mean, sd = sd)
})

data_facet <- as.data.frame(x = c("February", "March", "April","May","June","July","August","September"))
data_facet$inf <- info2017$mean - info2017$sd
data_facet$sup <- info2017$mean + info2017$sd
data_facet$xmin <- -Inf
data_facet$xmax <- Inf
colnames(data_facet) <- c("group","ymin","ymax","xmin","xmax")


# #NAVARRO FIG 4A - OLD
# b <- ggplot(perioddf2, aes(x = NORM_TOT_SIGMABIS, y = density, group = juld))+
#  geom_path(size = 0.5) + facet_grid(~group, scales = "free_x") +
#  scale_y_reverse() +   ylab("Potential density anomaly (kg/m³)") +
#  #xlab("Normalized chlorophyll a") +
#  geom_hline(data = data.frame(yint=as.vector(info2017$mean+info2017$sd),
#                               group=c("February", "March", "April","May","June","July","August","September")),
#                               aes(yintercept=yint), colour = "black", linetype="dashed") +
#  geom_hline(data = data.frame(yint=as.vector(info2017$mean-info2017$sd),
#                               group=c("February", "March", "April","May","June","July","August","September")),
#             aes(yintercept=yint), colour = "black",linetype="dashed") +
#  #geom_hline(yintercept = mean(info2017$mean), colour = "red") +
#  geom_rect(data=data_facet, aes(xmin=xmin, xmax=xmax,
#                                          ymin=ymin, ymax=ymax,
#                                          alpha=0,fill=factor("yellow")), inherit.aes = FALSE) +
#  scale_fill_manual(values = "yellow") + theme(legend.position="none") +
#  theme(text=element_text(size=12)) + #scale_x_continuous(breaks = c(0, 0.5))
#  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
#   axis.ticks.x=element_blank())

#NAVARRO FIG 4A
a <- ggplot(perioddf2, aes(x = NORM_TOT_SIGMABIS, y = density, group = juld))+
  geom_path(size = 0.5) + facet_grid(~group, scales = "free_x") +
  scale_y_reverse() +  ylab(expression(Potential~density~anomaly~(kg/m^3))) +
  #xlab("Normalized chlorophyll a") +
  geom_hline(yintercept = range_jan_feb2017[1], colour = "black", lty = "dashed") +
  geom_hline(yintercept = range_jan_feb2017[2], colour = "black", lty = "dashed") +
  annotate("rect", xmin=-Inf, xmax=+Inf, ymin = range_jan_feb2017[1], ymax = range_jan_feb2017[2],
           alpha = 0.5, fill="yellow", color = NA) +
  #scale_fill_manual(values = "yellow") + theme(legend.position="none") +
  theme(text=element_text(size=12)) + #scale_x_continuous(breaks = c(0, 0.5))
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


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
             DOY = tmp$DOY[1], juld = tmp$juld[1], month = tmp$month[1],
             year = tmp$year[1])
})

tmp <- TOTAL_DCM_PROFILES2[TOTAL_DCM_PROFILES2$year == YEARPAR,]

Per2 <- tmp[tmp$month == 2,]
Per3 <- tmp[tmp$month == 3,]
Per4 <- tmp[tmp$month == 4,]
Per5 <- tmp[tmp$month == 5,]
Per6 <- tmp[tmp$month == 6,]
Per7 <- tmp[tmp$month == 7,]
Per8 <- tmp[tmp$month == 8,]
Per9 <- tmp[tmp$month == 9,]

Per2$group <- c("February")
Per3$group <- c("March")
Per4$group <- c("April")
Per5$group <- c("May")
Per6$group <- c("June")
Per7$group <- c("July")
Per8$group <- c("August")
Per9$group <- c("September")

perioddf <- rbind(Per2, Per3, Per4, Per5, Per6, Per7, Per8, Per9)

#ordering facet
perioddf$group = factor(perioddf$group, 
                        levels=c("February", "March", "April","May","June","July","August","September"))

NORM_TOT_DCM_DEPTH <- normalize(perioddf$fluo)
perioddf2 <- cbind(perioddf, NORM_TOT_DCM_DEPTH)

perioddf2 <- perioddf2[!(perioddf2$month == 1 | perioddf2$month == 10 |
                           perioddf2$month == 11 | perioddf2$month == 12),]

rownames(perioddf2) <- NULL

#ADD MEAN DCM DEPTH AND PAR 1 & 10
dcm2017 <- ldply(as.list(2:9), function(i){
  tmp2 <- perioddf2[perioddf2$month == i,]
  tmp2 <- transform(tmp2,id=as.numeric(factor(juld)))
  dcm <- ddply(tmp2, ~juld, summarize,
               maxfluo = mean(depth[which.max(fluo)], na.rm = T))
  mean <- mean(dcm$maxfluo)
  sd <- sd(dcm$maxfluo)
  data.frame(mean = mean, sd = sd)
})

iso2017 <- ldply(as.list(2:9), function(i){
  tmp <- isolumes[isolumes$month == i,]
  p10 <- mean(tmp$iso10)
  p1 <- mean(tmp$iso1)
  data.frame(month = i, p10 = p10, p1 = p1)
})

#ADD MEAN MONTH MLD
mld2017 <- MLDdf[MLDdf$year == YEARPAR,]
mld2017 <- mld2017[!(mld2017$month == 1 | mld2017$month == 10 |
                       mld2017$month == 11 | mld2017$month == 12),]
mld2017 <- ddply(mld2017, ~month, summarize,
                 meanMLD = mean(MLD),
                 maxMLD = max(MLD))

#NAVARRO 4B
b <- ggplot(perioddf2, aes(x=NORM_TOT_DCM_DEPTH, y=depth, group=juld)) +
  geom_path(size = 0.5) + facet_grid(~group, scales = "free_x") +
  ylab("Depth (m)") +
  xlab("Normalized chlorophyll a") +
  theme(text=element_text(size=12)) + scale_x_continuous(breaks = c(0, 0.5)) + 
  scale_y_reverse(limits = c(80,0)) +
  geom_hline(data = data.frame(yint=as.vector(iso2017$p10), 
                               group=c("February", "March", "April","May","June","July","August","September")),
             aes(yintercept=yint), colour = "red") +
  geom_hline(data = data.frame(yint=as.vector(iso2017$p1), 
                               group=c("February", "March", "April","May","June","July","August","September")),
             aes(yintercept=yint), colour = "red", linetype = "dotted") +
  # geom_hline(data = data.frame(yint=as.vector(dcm2017$mean), 
  #                              group=c("February", "March", "April","May","June","July","August","September")),
  #            aes(yintercept=yint), colour = "yellow")  +
  geom_hline(data = data.frame(yint = as.vector(mld2017$meanMLD),
                               group = c("February", "March", "April","May","June","July","August","September")),
             aes(yintercept = yint), colour = "green") 

grid.arrange(a,b,nrow=2,ncol=1)

#NAVARRO 3A and 3B
tmp <- TOTAL_DCM_PROFILES2
NORM_TOT_DCM_DEPTH <- normalize(tmp$fluo)
tmp2 <- cbind(tmp, NORM_TOT_DCM_DEPTH)
tmp2 <- tmp2[!(tmp2$month == 1 | tmp2$month == 10 |
                 tmp2$month == 11 | tmp2$month == 12),]
rownames(tmp2) <- NULL


a <- ggplot(tmp2, aes(x = NORM_TOT_DCM_DEPTH, y = depth, group = juld, colour = month)) + 
  geom_path() + facet_grid(~year, scales = "free_x") + 
  ylab("Depth (m)") + xlab("Normalized chlorophyll a") + 
  theme(text=element_text(size=12)) + scale_x_continuous(breaks = c(0, 0.5)) + 
  scale_y_reverse(limits = c(80,0)) + 
  scale_color_gradientn(colors=rev(c("darkred", "red", "orange", "yellow",
                                     "green3", "lightgreen", "steelblue2", "steelblue4")),
                        breaks = c(2,3,4,5,6,7,8,9),
                        labels = c("February", "March", "April","May","June","July","August","September"),
                        guide = guide_colourbar(title = "", barheight = 10)) +
  theme(text=element_text(size=12)) + #scale_x_continuous(breaks = c(0, 0.5))
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank())



tmp <- TOTAL_SIGMA_PROFILES2
NORM_TOT_SIGMA <- normalize(tmp$fluo)
tmp2 <- cbind(tmp, NORM_TOT_SIGMA)
tmp2 <- tmp2[!(tmp2$month == 1 | tmp2$month == 10 |
                 tmp2$month == 11 | tmp2$month == 12),]
rownames(tmp2) <- NULL

rerangedf <- rangedf
rerangedf$xmin <- -Inf
rerangedf$xmax <- Inf
colnames(rerangedf) <- c("year","ymin","ymax","xmin","xmax")

b <- ggplot(tmp2, aes(x = NORM_TOT_SIGMA, y = density, group = juld, colour = month)) + 
  geom_path() + facet_grid(~year, scales = "free_x") + 
  ylab(expression(Potential~density~anomaly~(kg/m^3))) + xlab("Normalized chlorophyll a") + 
  theme(text=element_text(size=12)) + scale_x_continuous(breaks = c(0, 0.3)) + 
  scale_y_reverse() + 
  scale_color_gradientn(colors=rev(c("darkred", "red", "orange", "yellow",
                                     "green3", "lightgreen", "steelblue2", "steelblue4")),
                        breaks = c(2,3,4,5,6,7,8,9),
                        labels = c("February", "March", "April","May","June","July","August","September"),
                        guide = guide_colourbar(title = "", barheight = 10)) +
  # theme(text=element_text(size=13)) + #scale_x_continuous(breaks = c(0, 0.5)) 
  # theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  geom_hline(data=rangedf[rangedf$year == 2014,], aes(yintercept=rangedf$inf[1]), colour="black", lty = "dashed") +
  geom_hline(data=rangedf[rangedf$year == 2015,], aes(yintercept=rangedf$inf[2]), colour="black", lty = "dashed") + 
  geom_hline(data=rangedf[rangedf$year == 2016,], aes(yintercept=rangedf$inf[3]), colour="black", lty = "dashed") +
  geom_hline(data=rangedf[rangedf$year == 2017,], aes(yintercept=rangedf$inf[4]), colour="black", lty = "dashed") +
  geom_hline(data=rangedf[rangedf$year == 2014,], aes(yintercept=rangedf$sup[1]), colour="black", lty = "dashed") +
  geom_hline(data=rangedf[rangedf$year == 2015,], aes(yintercept=rangedf$sup[2]), colour="black", lty = "dashed") +
  geom_hline(data=rangedf[rangedf$year == 2016,], aes(yintercept=rangedf$sup[3]), colour="black", lty = "dashed") +
  geom_hline(data=rangedf[rangedf$year == 2017,], aes(yintercept=rangedf$sup[4]), colour="black", lty = "dashed") +
  geom_rect(data=rerangedf, aes(xmin=xmin, xmax=xmax,
                                ymin=ymin, ymax=ymax),
            alpha=.2,fill=factor("black"), inherit.aes = FALSE)

grid.arrange(a,b, ncol = 1, nrow = 2)

######### DCM avec la PAR (LETELIER 2004)

PARanalysis <-ldply(as.list(PAR_files),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = FALSE, suppress_dimvals = FALSE)
  
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
  
  show(file)
  
  chladf <- ExtractVar("CHLA", FloatInfo)
  chladf <- chladf[-(which(chladf$qc == 4)),]
  pardf <- ExtractVar("DOWNWELLING_PAR", FloatInfo)
  tempdf <- ExtractVar("TEMP",FloatInfo)
  psaldf <- ExtractVar("PSAL",FloatInfo)
  
  chladf <- var_average(chladf)
  pardf <- var_average(pardf)
  tempdf <- var_average(tempdf)
  tempdf$value[which(tempdf$qc == 4 | tempdf$qc == 3)] <- NA
  psaldf <- var_average(psaldf)
  psaldf$value[which(psaldf$qc == 4 | psaldf$qc == 3)] <- NA
  pardf <- match_df(pardf, chladf, on = c("juld","depth"))
  chladf <- match_df(chladf, pardf, on = c("juld","depth"))
  tempdf <- match_df(tempdf, chladf, on = c("juld","depth"))
  psaldf <- match_df(psaldf, pardf, on = c("juld","depth"))
  chladf <- match_df(chladf, psaldf, on = c("juld","depth"))
  pardf <- match_df(pardf, psaldf, on = c("juld","depth"))
  
  psal <- gsw_SA_from_SP(psaldf$value,psaldf$depth,psaldf$lon,psaldf$lat)
  temp <- gsw_CT_from_t(psal,tempdf$value,psaldf$depth)
  sigma <- gsw_sigma0(psal,temp)
  
  #Construction du data frame final a lieu ici
  data.frame(depth         = chladf$depth,
             juld          = chladf$juld,
             fluo          = chladf$value,
             par           = pardf$value,
             sigma         = sigma,
             day           = month.day.year(chladf$juld,c(1,1,1950))$day,
             month         = month.day.year(chladf$juld,c(1,1,1950))$month,
             year          = month.day.year(chladf$juld,c(1,1,1950))$year,
             DOY           = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon           = chladf$lon,
             lat           = chladf$lat,
             Platform      = as.numeric(unique(id)),
             type          = "Argo")
})

PARanalysis <- transform(PARanalysis,id=as.numeric(factor(juld)))

tmp <- PARanalysis[PARanalysis$id == i,]
ggplot(tmp, aes(x = fluo, y = par)) + 
  geom_path() + ggtitle(i)
i <- i + 1 

parreject <- c(9,10,12,16,26,44,47,55,72,83,101,45,43,42,40,39,37,35,33,32,48,61,85,
               208,206,192,152,148,141,126,123,110,116,137)

PARanalysis <- PARanalysis[!(PARanalysis$id %in% parreject),]

tmp <- PARanalysis[PARanalysis$year == 2017,]
tmp <- PARanalysis
tmp <- tmp[-(which(is.na(tmp$sigma)==TRUE)),]
rownames(tmp) <- NULL
tmp <- transform(tmp,id=as.numeric(factor(juld)))

df <- ddply(tmp, ~juld, summarize,
            bottom = max(depth, na.rm =T),
            min = min(depth, na.rm = T),
            sigmamax = sigma[which.max(fluo)],
            par = par[which.max(fluo)],
            maxfluo = max(fluo))

df <- df[-(which(df$bottom < 100)),]
df <- df[-(which(df$bottom == df$min)),]
df <- df[-(which(df$min > 10)),]

tmp <- tmp[tmp$juld %in% df$juld,]

range2017 <- range(info2017$mean)

ggplot(tmp, aes(x = fluo, y = par)) +
  geom_point(size = 1) + geom_point(data = df, aes(x = maxfluo, y = par), color = "red", size = 1) +
  xlab(expression(Chlorophyll~a~(mg/m^3))) + #ylab("PAR (microMoleQuanta/m²/sec)")
  #ylab(expression(mu),expression(mol~quanta~(m^-2~s^-1))) +  
  theme(text=element_text(size=12)) +
  ylab (expression(paste(
    "PAR (",
    mu, mol," ","quanta ",m^-2, s^-1,
    ")", sep="")))

ggplot(tmp, aes(x = fluo, y = par)) +
  geom_point(size = 1) + geom_point(data = df, aes(x = maxfluo, y = par), color = "red", size = 1) +
  xlab(expression(Chlorophyll~a~(mg/m^3))) + #ylab("PAR (microMoleQuanta/m²/sec)")
  #ylab(expression(mu),expression(mol~quanta~(m^-2~s^-1))) +  
  theme(text=element_text(size=12)) +
  ylab (expression(paste(
    "PAR (",
    mu, mol," ","quanta ",m^-2, s^-1,
    ")", sep=""))) + ylim(c(0,300)) + geom_hline(yintercept = 5.787, lty = "dashed") + xlim(c(0,2.5))

#PRE-DIVA
#get MOD DCM 
MODDCM <- TOTAL2[TOTAL2$Category == "MODDCM",]

a <- ggplot(TOTAL2, aes(x= DOY, y = sigma_dcm)) +
  geom_point() + scale_color_manual(values = c("black", "red")) + theme(legend.position = "none") +
  geom_smooth() + ylab(expression(DCM~potential~density~anomaly~(kg/m^3))) + xlab("Day Of Year (DOY)") + 
  geom_point(data = MODDCM, aes(x = DOY, y = sigma_dcm), colour = "red") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(text=element_text(size=12)) 

b <- ggplot(TOTAL2, aes(x= DOY, y = Zmax)) +
  geom_point() + scale_color_manual(values = c("black", "red")) + theme(legend.position = "none") +
  geom_smooth() + ylab("DCM depth (m)") + xlab("Day Of Year (DOY)") + 
  geom_point(data = MODDCM, aes(x = DOY, y = Zmax), colour = "red") + theme(text=element_text(size=12)) 

grid.arrange(a,b,nrow = 2, ncol = 1)

###########

#mean MLD from 2014 to 2017

check3 <- ddply(DCM, ~year~month, summarize,
                n = length(month),
                meandepth = mean(Zmax),
                juld = juld[1])

new_time <- as.Date(check3$juld, origin = as.Date("1950-01-01"))
new_time <- format(as.Date(new_time), "%Y-%m")
check3$date <- new_time
#check3$date <- as.POSIXct(check3$date)
#library(scales)

# ggplot(check3, aes(date)) +
#   geom_bar(aes(weight=n)) + ylab("DCM mean depth (m)") + xlab("Time") + 
#   theme(text=element_text(size=12)) + theme(axis.text.x = element_text(color="black", size=10, angle=45)) +
#   scale_x_datetime(date_breaks = "3 months")


ggplot(check3, aes(date)) +
  geom_bar(aes(weight=meandepth), fill="steelblue4", colour="black") + ylab("DCM mean depth (m)") + xlab("Time") +
  theme(text=element_text(size=12)) + theme(axis.text.x = element_text(color="black", size=10, angle=90)) +
  geom_bar(data = check3, aes(weight=n, colour="black"), fill="steelblue3",colour="black")


