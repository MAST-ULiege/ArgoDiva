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

# Extracting fluorescence values (CHLA and CDOM) from metadata files
# meta_filename <- c("6900807_meta.nc", "6901866_meta.nc", "7900591_meta.nc", "7900592_meta.nc")
# Extract_PREDEPLOYMENT_CALIB_COEFFICIENT <- ldply(as.list(meta_filename), function(file){
#   ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
#   test <- ncvar_get(ncfile,"PREDEPLOYMENT_CALIB_COEFFICIENT")
# })

param_6900807 <- c(0.0072, 47) #SCALE_CHLA, DARK_CHLA
param_6901866 <- c(0.0073, 49)
param_7900591 <- c(0.0121, 48)
param_7900592 <- c(0.0121, 50)
fluo_param <- list(param_6900807, param_6901866, param_7900591,
                   param_7900592)

# SECTION 1 : Extracting data from CHLA-equipped ARGO floats ------------------
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

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

#Get each profile from ARGO data -----
profiledf <- ldply(as.list(filename),function(file){
  
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
  
  chladf <- ExtractVar("CHLA",FloatInfo)
  
  # CHLA = (FLUO_CHLA - DARK_CHLA) * SCALE_CHLA
  # => FLUO_CHLA = (CHLA / SCALE_CHLA) + DARK_CHLA
  
  #fluorescence conversion, respectively to each float
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
  }
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
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

#Save
profiledf_raw <- profiledf
profiledf_raw <- transform(profiledf_raw,id=as.numeric(factor(juld)))

#Remove bad data 
bad_data <- which(profiledf$qc == 4)
profiledf$fluo[bad_data] <- NA

#Test median filter
#mav <- function(x, n=7){filter(x, rep(1/n,n), sides=2)}

#Set profiles id
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

#Construction of an ARGO-only dataframe -> first guess value for Zmax [NLS fit]
argodf<-ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
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
  
  chladf <- ExtractVar("CHLA", FloatInfo)
  
  #fluorescence conversion, respectively to each float
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
  }
  
  #Chla subset df et final df
  subchladf <- subset(chladf,select=c("depth","juld","value","qc","lon","lat"))
  colnames(subchladf)[which(colnames(subchladf)=="value")]<-"FLUO"
  
  chladf <- ddply(subchladf,~juld,summarize,
                  qc = qc[which.max(FLUO)],
                  depthmax = depth[which.max(FLUO)],
                  maxvalue = FLUO[which.max(FLUO)],
                  integration = sum(FLUO),
                  bottomdepth = max(depth),
                  lon=mean(lon),
                  lat=mean(lat))#mean car la bouée a bougé sur la journée ?
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
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

#Assign a profile number according to the (unique due to two decimals) julian day 
argodf <- transform(argodf,id=as.numeric(factor(juld)))

#Get profiles id that do not answer to some criteria
crit1 <- argodf[(argodf$bottom < 100),]#Remove profiles too shallow -> 2013 anyway et celui de 2015 sera viré sar stuck
crit2 <- argodf[(argodf$depthmax >= 100),]#NLS is performed on the first 100m (max > 100m should not happen anyway)
crit3 <- argodf[(argodf$depthmax == argodf$bottom),]#'Stuck' profiles (all values at the same depth), NOTE: Ils ne sont pas à QC = 4? Weird
criteriadf <- rbind(crit1, crit2, crit3)
criteriadf <- unique(criteriadf)

#Save
argodf_raw <- argodf

#remove profiles that have been categorized as bad profiles
for (i in 1:length(criteriadf$id)){
  argodf <- argodf[!(argodf$id == criteriadf$id[i]),]
}

argodf <- transform(argodf,id=as.numeric(factor(juld)))

#REMOVE 2013 ET 2018 (because not FULL year)
argodf <- argodf[!argodf$year == 2013,]
argodf <- argodf[!argodf$year == 2018,]
#REMOVE 2014 ? -> Not enough data?

for (i in 1:length(criteriadf$id)){
  profiledf <- profiledf[!(profiledf$id == criteriadf$id[i]),]
}

#remove profiles from 2013 and 2018
profiledf <- profiledf[!profiledf$year == 2013,]
profiledf <- profiledf[!profiledf$year == 2018,]

##new id before gaussian elimination
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

# MLD 

densitydf<-ldply(as.list(filename),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
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
  
  tempdf <- ExtractVar("TEMP_ADJUSTED",FloatInfo)
  if (all(is.na(tempdf)) == TRUE) {
    tempdf <- ExtractVar("TEMP",FloatInfo)
  }
  
  bad_data <- which(tempdf$qc == 4)
  tempdf$value[bad_data] <- NA
  bad_data <- which(tempdf$qc == 3)
  tempdf$value[bad_data] <- NA
  
  psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    psaldf <- ExtractVar("PSAL",FloatInfo)
  }
  
  bad_data <- which(psaldf$qc == 4)
  psaldf$value[bad_data] <- NA
  bad_data <- which(psaldf$qc == 3)
  psaldf$value[bad_data] <- NA
  
  psal <- gsw_SA_from_SP(psaldf$value,psaldf$depth,psaldf$lon,psaldf$lat)
  temp <- gsw_CT_from_t(psal,tempdf$value,tempdf$depth)
  rho_anomaly <- gsw_sigma0(psal,temp)
  
  data.frame(density_anomaly = rho_anomaly,
             depth = tempdf$depth, juld = tempdf$juld, lon = tempdf$lon,
             lat = tempdf$lon,
             day             = month.day.year(tempdf$juld,c(1,1,1950))$day,
             month           = month.day.year(tempdf$juld,c(1,1,1950))$month,
             year            = month.day.year(tempdf$juld,c(1,1,1950))$year,
             DOY             = as.integer(strftime(as.Date(tempdf$juld,origin = '1950-01-01'), 
                                                   format ="%j"))#Day Of Year)
  )
  
})

densitydf <- transform(densitydf,id=as.numeric(factor(juld)))

for (i in 1:length(criteriadf$id)){
  densitydf <- densitydf[!(densitydf$id == criteriadf$id[i]),]
}

#remove profiles from 2013 and 2018
densitydf <- densitydf[!densitydf$year == 2013,]
densitydf <- densitydf[!densitydf$year == 2018,]

densitydf <- transform(densitydf,id=as.numeric(factor(juld)))

MLDdf <- ldply(as.list(1:length(unique(densitydf$id))), function(i){#le compiler plusieurs fois à la main puis le global fonctionne...... WHY??????
# MLDdf <- ldply(as.list(750:839), function(i){
    tmp <- densitydf[densitydf$id == i,]
  
  if(length(tmp$density_anomaly[!is.na(tmp$density_anomaly)==TRUE])>=2) {
    if (min(tmp$depth > 10)){
      rho_surf <- tmp$density_anomaly[which(tmp$density_anomaly == min(tmp$density_anomaly, na.rm = T))]
      MLD <- max(tmp$depth[tmp$density_anomaly <= (rho_surf + 0.03)],na.rm = T)
    } else{
      rho_surf <- NA
      rho_surf <- approx(tmp$depth, tmp$density_anomaly, 10)$y
    }
    
    if (is.na(rho_surf)==FALSE) {
      MLD <- max(tmp$depth[tmp$density_anomaly <= (rho_surf + 0.03)],na.rm = T)
    }
  }else if(length(tmp$density_anomaly[!is.na(tmp$density_anomaly)==TRUE])==0){
    rho_surf <- NA
    MLD <- NA
  }else{
    rho_surf <- NA
    MLD <- NA
  }
  
  data.frame(rho_surf = rho_surf, MLD = MLD, juld = tmp$juld[1], day = month.day.year(tmp$juld[1],c(1,1,1950))$day,
             month = month.day.year(tmp$juld[1],c(1,1,1950))$month,
             year = month.day.year(tmp$juld[1],c(1,1,1950))$year,
             DOY = as.integer(strftime(as.Date(tmp$juld[1],origin = '1950-01-01'), 
                                       format ="%j")))#Day Of Year
  # i <- i + 1

  })

saveMLDdf <- MLDdf
MLDdf <- transform(MLDdf,id=as.numeric(factor(juld)))

#remove profiles for which rho_surf ou MLD == NA
# TO DO? -> au pire, pas de correction de quenching -> FAUX DCM?
remove_uncertain_profiles <- which(is.na(MLDdf$rho_surf))
for (i in 1:length(remove_uncertain_profiles)){
profiledf <- profiledf[!(profiledf$id == remove_uncertain_profiles[i]),]
MLDdf <- MLDdf[!(MLDdf$id == remove_uncertain_profiles[i]),]
}

profiledf <- transform(profiledf,id=as.numeric(factor(juld)))
MLDdf <- transform(MLDdf,id=as.numeric(factor(juld)))

NPQcorrection<- function(MLD, tmp) {
  
  f <- tmp$fluo[!is.na(tmp$fluo) & tmp$depth <= MLD]
  d <- tmp$depth[!is.na(tmp$fluo) & tmp$depth <= MLD]
  zMax <- d[which.max(f)]
  Max <- max(f)
  Corfluo <- tmp$fluo
  
  if(min(f[d<zMax])<(0.9*Max)) {#find the justification for this -> cela veut dire qu'il n'y a probablement pas de NPQ
    Corfluo[tmp$depth<zMax] <- Max
  }
  
  return(Corfluo)
}

fluo_adjusted <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id == i,]
  MLD <- MLDdf$MLD[i]
  fluo_adjusted <- NPQcorrection(MLD, tmp)
  
  data.frame(fluo_adjusted = fluo_adjusted, id = i)  
})

profiledf <- profiledf[order(profiledf$id),] 
profiledf <- cbind(profiledf, fluo_adjusted)


i <-65
tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex])
x <- tab$x
tab2 <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo_adjusted[1:depthindex])
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab2$x, tab2$y,col=2, add=T, xlim=range(tab2$x), lwd = 2)
i <- i + 1


#Graphe reconstruction
# i <- 166
# tmp <- profiledf[profiledf$id==i,]
# depthindex <- which.min(tmp$depth <= 100)
# tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
# #tab <- data.frame(x=tmp$depth,y=tmp$chla)
# plot(y~x, data=tab, type="l", main = i)
# tmp2 <- approx(tab$x,tab$y, method="linear")
# #abline(v=seq(0,900,100) , col="grey" , lwd=0.6)
# lines(y~x, data=tmp2, type="l", col="red")
# i <- i+1


#HERE PROFILES SHOULD BE OK BECAUSE BAD DATA HAVE BEEN REMOVED#

#### Mignot 2011 ####

fgauss <- function(x, Fsurf, Zdemi, Fmax, Zmax, dz){
  Fsurf*exp((-log(2)/Zdemi)*x) + Fmax*exp(-(x-Zmax)^2/dz^2)
}

fsigmoid <- function(x, Fsurf, Zdemi, s){
  Fsurf*(1/(1+exp((Zdemi-x)*s)))
}

gaussiandf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  #gaussiandf <- ldply(as.list(176:185), function(i){
    
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  
  depthindex <- which.min(tmp$depth <= 100)
  off_fluo <- tmp$fluo_adjusted[which.min(tmp$fluo_adjusted[1:depthindex])]
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo_adjusted[1:depthindex]-off_fluo)
  maxindex <- which.max(tmp$fluo_adjusted)
  #Parameters estimation
  nonNAindex <- which(!is.na(tmp$fluo_adjusted))
  firstnonNA <- nonNAindex[1]
  Fsurf <- tmp$fluo_adjusted[firstnonNA]
  Zmax <- tmp$depth[maxindex]
  if(Zmax > 100){
    Zmax <- 50
  }
  Fmax <- tmp$fluo_adjusted[maxindex]
  Zdemi <- tmp$depth[which(tmp$fluo_adjusted <= Fsurf/2)[2]]
  if(is.na(Zdemi)){#code à la bourrain
    indexZ <- which(tmp$fluo_adjusted < Fsurf)
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

sigmoidf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  depthindex <- which.min(tmp$depth <= 100)
  off_fluo <- tmp$fluo_adjusted[which.min(tmp$fluo_adjusted[1:depthindex])]
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo_adjusted[1:depthindex]-off_fluo)
  maxindex <- which.max(tmp$fluo_adjusted)
  #Parameters estimation
  nonNAindex <- which(!is.na(tmp$fluo_adjusted))
  firstnonNA <- nonNAindex[1]
  Fsurf <- tmp$fluo_adjusted[firstnonNA]
  Zdemi <- tmp$depth[which(tmp$fluo_adjusted <= Fsurf/2)[1]]
  if(is.na(Zdemi)){
    indexZ <- which(tmp$fluo_adjusted < Fsurf)
    Zdemi <- tmp$depth[indexZ[which.max(indexZ > maxindex)]]
  }
  i1 <- which.min(tmp$depth <= Zdemi - 5)
  i2 <- which.min(tmp$depth <= Zdemi + 5)
  d1 <- tmp$depth[i1]
  d2 <- tmp$depth[i2]
  s <- -(Fsurf/2)/(d2-d1)
  
  if (!(Fsurf > 3*tmp$fluo_adjusted[length(tmp$fluo_adjusted)])){
    s <- -2.5
  }
  
  
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

#visualisation triplot
i <- 2
tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
off_fluo <- tmp$fluo_adjusted[which.min(tmp$fluo_adjusted[1:depthindex])]
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo_adjusted[1:depthindex]-off_fluo)
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
                          sigmoidf$s[i]),col=2, add=T, xlim=range(tab$x), lwd = 2)
lines(x = tab$x, fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
                        gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i + 1

#Coefficient of determination
Rcoefdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  depthindex <- which.min(tmp$depth <= 100)
  off_fluo <- tmp$fluo_adjusted[which.min(tmp$fluo_adjusted[1:depthindex])]
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo_adjusted[1:depthindex]-off_fluo)
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

#classification des profiles

R_crit <- 0.8 #80% de variance expliquée (Navarro c'était 90%)

classifdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  
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
ngauss <- length(which(classifdf$classif == "gaussian"))/length(unique(profiledf$id))*100
nsigmoid <- length(which(classifdf$classif == "sigmoid"))/length(unique(profiledf$id))*100
nother <- length(which(classifdf$classif == "other"))/length(unique(profiledf$id))*100

#gaussian profile (classified) before a new criteria elimination
gaussprofiles <- classifdf[classifdf$classif == "gaussian",]

#remove profiles where depthmax (DCM) < 1m ===> navarro
gaussprofiles <- gaussprofiles[gaussprofiles$depthmax > 1,]

#Navarro critère 1 is back  ! Remove DCM (from gaussian designated profiles) below 1m
#actually, I think it would be interesting to remove DCM peut-être encore < 5m..
#purement subjectif

gaussprofiles <- transform(gaussprofiles,id=as.numeric(factor(juld)))

#NOW FLUO IS REFERRED TO ADJUSTED FLUO (PART 1 IS DONE)

#get pseudo gaussian profiles
gaussianprofdf <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- profiledf[profiledf$juld == gaussprofiles$juld[i],]
  data.frame(depth = tmp$depth, juld = tmp$juld, fluo = tmp$fluo_adjusted, qc = tmp$qc,
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


#VISU analysis of Modified_DCM
#10% of profiles were considered as non modified_DCM

j <- 1
i <- t4[j]
tmp <- profiledf[profiledf$juld == gaussprofiles$juld[i],]
depthindex <- which.min(tmp$depth <= 100)
off_fluo <- tmp$fluo_adjusted[which.min(tmp$fluo_adjusted[1:depthindex])]
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo_adjusted[1:depthindex])
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
  #theme(axis.text.x = element_text(angle=0, vjust=2))

modified_DCM <- temporal_stats[temporal_stats$Category == "Modified_DCM",]
modified_DCM <- transform(modified_DCM,id=as.numeric(factor(juld)))
DCM <- temporal_stats[temporal_stats$Category == "DCM",]
DCM <- transform(DCM,id=as.numeric(factor(juld)))

#get DCM profiles
DCMprofiles <- ldply(as.list(1:length(DCM$juld)), function(i){
  tmp <- gaussprofiles[gaussprofiles$juld == DCM$juld[i],]
  tmp3 <- gaussdata[gaussdata$juld == DCM$juld[i],]
  tmp2 <- profiledf[profiledf$juld == DCM$juld[i],]
  n <- length(tmp2$depth)
  data.frame(depth = tmp2$depth, juld = tmp2$juld, fluo = tmp2$fluo_adjusted, qc = tmp2$qc,
             day = tmp2$day, month = tmp2$month, year = tmp2$year, DOY = tmp2$DOY,
             Platform = tmp2$Platform, lon = tmp2$lon, lat = tmp2$lat, 
             Fsurf = rep(tmp3$Fsurf, n), Fmax = rep(tmp3$Fmax, n), Zmax = rep(tmp3$Zmax, n),
             dz = rep(tmp3$dz, n), Zdemi = rep(tmp3$Zdemi, n))
})

#get modified DCM
modified_DCMprofiles <- ldply(as.list(1:length(modified_DCM$juld)), function(i){
  tmp <- gaussprofiles[gaussprofiles$juld == modified_DCM$juld[i],]
  tmp3 <- gaussdata[gaussdata$juld == modified_DCM$juld[i],]
  tmp2 <- profiledf[profiledf$juld == modified_DCM$juld[i],]
  n <- length(tmp2$depth)
  data.frame(depth = tmp2$depth, juld = tmp2$juld, fluo = tmp2$fluo_adjusted, qc = tmp2$qc,
             day = tmp2$day, month = tmp2$month, year = tmp2$year, DOY = tmp2$DOY,
             Platform = tmp2$Platform, lon = tmp2$lon, lat = tmp2$lat, 
             Fsurf = rep(tmp3$Fsurf, n), Fmax = rep(tmp3$Fmax, n), Zmax = rep(tmp3$Zmax, n),
             dz = rep(tmp3$dz, n), Zdemi = rep(tmp3$Zdemi, n))
})

DCMprofiles <- transform(DCMprofiles,id=as.numeric(factor(juld)))
modified_DCMprofiles <- transform(modified_DCMprofiles, id=as.numeric(factor(juld)))

#Check Visu DCM
tmp <- DCMprofiles[DCMprofiles$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
off_fluo <- tmp$fluo[which.min(tmp$fluo[1:depthindex])]
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex]-off_fluo)
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, fgauss(x,tmp$Fsurf[i], tmp$Zdemi[i],
                        tmp$Fmax[i], tmp$Zmax[i], tmp$dz[i]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i + 1

#Check Visu modified_DCM
tmp <- modified_DCMprofiles[modified_DCMprofiles$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
off_fluo <- tmp$fluo[which.min(tmp$fluo[1:depthindex])]
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex]-off_fluo)
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, fgauss(x,tmp$Fsurf[i], tmp$Zdemi[i],
                        tmp$Fmax[i], tmp$Zmax[i], tmp$dz[i]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i + 1

#fusion density et depth R codes -> need the find the density associated to the depth DCM....

#function (code pourri -> bcp de problèmes bizarres.......... /-|-\" 

#reorder and sort by filename for rho_DCM
DCM <- DCM[order(DCM$file),]
modified_DCM <- modified_DCM[order(modified_DCM$file),]

#datadf is either DCM or modified_DCM
rho_dcm_from_depth_dcm <- function(file, borne_inf, borne_sup, datadf){
  
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
  
  tempdf <- ExtractVar("TEMP_ADJUSTED",FloatInfo)
  if (all(is.na(tempdf)) == TRUE) {
    tempdf <- ExtractVar("TEMP",FloatInfo)
  }
  
  bad_data <- which(tempdf$qc == 4)
  tempdf$value[bad_data] <- NA
  bad_data <- which(tempdf$qc == 3)
  tempdf$value[bad_data] <- NA
  
  psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    psaldf <- ExtractVar("PSAL",FloatInfo)
  }
  
  bad_data <- which(psaldf$qc == 4)
  psaldf$value[bad_data] <- NA
  bad_data <- which(psaldf$qc == 3)
  psaldf$value[bad_data] <- NA
  
  df <- ldply(as.list(borne_inf:borne_sup), function(i){
    
    depth_dcm <- datadf$Zmax[i]
    julian_day <- datadf$juld[i]
    
    psaltmp <- psaldf[psaldf$juld == julian_day,]
    temptmp<- tempdf[tempdf$juld == julian_day,]
    
    indexdepth <- which.min(abs(temptmp$depth - depth_dcm))
    psaltmp <- psaltmp[indexdepth,]
    temptmp <- temptmp[indexdepth,]
    
    psal <- gsw_SA_from_SP(psaltmp$value,depth_dcm,psaltmp$lon,psaltmp$lat)
    temp <- gsw_CT_from_t(psal,temptmp$value,depth_dcm) 
    rho_anomaly <- gsw_sigma0(psal,temp)
    #i <- i + 1
    data.frame(rho_dcm = rho_anomaly, filename = datadf$file[i],
               lon = datadf$lon[i], lat = datadf$lat[i],
               day = datadf$day[i], month = datadf$month[i], year = datadf$year[i],
               DOY = datadf$DOY[i], juld = datadf$juld[i])
  })
}

gaussdata_raw <- gaussdata

#Choose your data
datadf <- modified_DCM
datadf <- DCM

index_dcm <- which(duplicated(datadf$file) == FALSE)

rho_dcm1 <- rho_dcm_from_depth_dcm(datadf$file[index_dcm[1]], index_dcm[1], index_dcm[2]-1, datadf)
rho_dcm2 <- rho_dcm_from_depth_dcm(datadf$file[index_dcm[2]], index_dcm[2], index_dcm[3]-1, datadf)
rho_dcm3 <- rho_dcm_from_depth_dcm(datadf$file[index_dcm[3]], index_dcm[3], index_dcm[4]-1, datadf)
rho_dcm4 <- rho_dcm_from_depth_dcm(datadf$file[index_dcm[4]], index_dcm[4], length(datadf$juld), datadf)
rho_dcm <- rbind(rho_dcm1, rho_dcm2, rho_dcm3, rho_dcm4)

rhoNA <- which(is.na(rho_dcm$rho_dcm) == TRUE)
modified_dcm_juld_to_remove <- rho_dcm$juld[rhoNA]#for below
dcm_juld_to_remove <- rho_dcm$juld[rhoNA]#for below
#remove NA's
rho_dcm <- rho_dcm[complete.cases(rho_dcm),]
rho_dcm <- transform(rho_dcm,id=as.numeric(factor(juld)))

rhoDCM <- rho_dcm
rho_modified_DCM <- rho_dcm 

# FRACTION OF VARIANCE UNEXPLAINED
# https://en.wikipedia.org/wiki/Fraction_of_variance_unexplained

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

#Fonction qui renvoie un data frame basé sur l'anomalie de densité verticale (y axis)
ExtractDensity <- function(chladf, FloatInfo){
  
  #Extraction of temperature and salinity profiles in order to
  #retrieve the potential density anomaly referenced to 0dbar  
  
  #Retrieve temp and psal in front of chlorophyll-a (
  #also see : ncvar_get(ncfile,"STATION_PARAMETERS"))
  #ncvar_get(ncfile,"VERTICAL_SAMPLING_SCHEME")
  
  tempdf <- ExtractVar("TEMP_ADJUSTED",FloatInfo)
  if (all(is.na(tempdf)) == TRUE) {
    tempdf <- ExtractVar("TEMP",FloatInfo)
  }
  
  tempdf <- var_average(tempdf)#for averaging depth duplicates for temperature
  #instead of removing bad lines, put NA's ---> will impact graphics
  bad_data <- which(tempdf$qc == 4)
  tempdf$value[bad_data] <- NA
  bad_data <- which(tempdf$qc == 3)
  tempdf$value[bad_data] <- NA
  
  psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    psaldf <- ExtractVar("PSAL",FloatInfo)
  }
  psaldf <- var_average(psaldf)#for averaging depth duplicates for salinity
  bad_data <- which(psaldf$qc == 4)
  psaldf$value[bad_data] <- NA
  bad_data <- which(psaldf$qc == 3)
  psaldf$value[bad_data] <- NA
  
  # MEME PLUS BESOIN DE ÇA !
  # #Extract matching rows temperature and salinity data frames
  # tempdf <- match_df(tempdf, psaldf, on = c("depth", "juld", "aprofile", "alevel"))
  # psaldf <- match_df(psaldf, tempdf, on = c("depth", "juld", "aprofile", "alevel"))
  # tempdf <- match_df(tempdf, chladf, on = c("depth", "juld"))
  # psaldf <- match_df(psaldf, chladf, on = c("depth", "juld"))
  # chladf <- match_df(chladf, psaldf, on = c("depth", "juld"))#135 lignes de plus que ce à quoi je m'attends...
  # #24101 pour tempdf et psaldf vs 24236 pour chladf et finaldf (voir après)
  # 
  #Sub data frames

  subtempdf2 <- subset(tempdf, select = c("value","depth","juld","lon","lat"))
  colnames(subtempdf2)[which(colnames(subtempdf2)=="value")]<-"TEMP"
  subpsaldf2 <- subset(psaldf, select = c("value","depth","juld"))
  colnames(subpsaldf2)[which(colnames(subpsaldf2)=="value")]<-"PSAL"
  #subchladf2 <- subset(chladf, select = c("fluo_adjusted","depth","juld","qc"))
  subchladf2 <- subset(chladf, select = c("fluo","depth","juld","qc"))
  colnames(subchladf2)[which(colnames(subchladf2)=="fluo")]<-"FLUO_adjusted"
  joindf <- join(subtempdf2,subpsaldf2, by = c("depth","juld"))
  subchladf2 <- match_df(subchladf2,joindf, on = c("depth", "juld"))
  finaldf <- join(subchladf2,joindf, by = c("depth","juld"))
  
  #use of gsw package
  #NEED CONVERSION OF PRACTICAL SALINITY TO ABSOLUTE SALINITY BEFORE USING THE FOLLOWING FUNCTION
  #NEED CONVERSION OF IN-SITU TEMPERATURE TO CONSERVATIVE TEMPERATURE
  psal <- gsw_SA_from_SP(finaldf$PSAL,finaldf$depth,finaldf$lon,finaldf$lat)
  temp <- gsw_CT_from_t(psal,finaldf$TEMP,finaldf$depth)
  rho_anomaly <- gsw_sigma0(psal,temp)
  
  rhodf <- data.frame(depth         = finaldf$depth,
                      rho_anomaly   = rho_anomaly,  
                      FLUO_ADJUSTED = finaldf$FLUO_adjusted,
                      TEMP          = temp,
                      PSAL          = psal,
                      juld          = finaldf$juld,
                      day           = month.day.year(finaldf$juld,c(1,1,1950))$day,
                      month         = month.day.year(finaldf$juld,c(1,1,1950))$month,
                      year          = month.day.year(finaldf$juld,c(1,1,1950))$year,
                      lon           = finaldf$lon,
                      lat           = finaldf$lat)
}

rho_DCM_modified_densityprofiledf<-ldply(as.list(filename),function(file){
  
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
             density       = densitydf$rho_anomaly,
             juld          = densitydf$juld,
             FLUO_ADJUSTED = densitydf$FLUO_ADJUSTED,
             day           = month.day.year(densitydf$juld,c(1,1,1950))$day,
             month         = month.day.year(densitydf$juld,c(1,1,1950))$month,
             year          = month.day.year(densitydf$juld,c(1,1,1950))$year,
             DOY           = as.integer(strftime(as.Date(densitydf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon           = densitydf$lon,
             lat           = densitydf$lat,
             Platform      = as.numeric(unique(id)),
             type          = "Argo")
})

rho_DCM_densityprofiledf_raw <- rho_DCM_densityprofiledf
rho_DCM_modified_densityprofiledf_raw <- rho_DCM_modified_densityprofiledf

#pour supprimer les doublons qui apparaissent
rho_DCM_densityprofiledf <- unique(rho_DCM_densityprofiledf)
rho_DCM_modified_densityprofiledf <- unique(rho_DCM_modified_densityprofiledf)
#Assign a profile number according to the (unique due to two decimals) julian day 
rho_DCM_densityprofiledf <- transform(rho_DCM_densityprofiledf,id=as.numeric(factor(juld)))
rho_DCM_modified_densityprofiledf <- transform(rho_DCM_modified_densityprofiledf,id=as.numeric(factor(juld)))

#REMOVE PROFILES FOR WHICH JULD HAS TO BE REMOVED BECAUSE RHO COULD NOT 
# BE COMPUTED WITH CERTAINTY
# CHECK THAT AFTER IT IS WELL 228 AND 183

#DCM cleaning
for (i in 1:length(dcm_juld_to_remove)){
  rho_DCM_densityprofiledf <- rho_DCM_densityprofiledf[!(rho_DCM_densityprofiledf$juld == dcm_juld_to_remove[i]),]
  DCMprofiles <- DCMprofiles[!(DCMprofiles$juld == dcm_juld_to_remove[i]),]
  DCM <- DCM[!DCM$juld == dcm_juld_to_remove[i],]
  }

rho_DCM_densityprofiledf <- transform(rho_DCM_densityprofiledf,id=as.numeric(factor(juld)))
DCMprofiles <- transform(DCMprofiles, id = as.numeric(factor(juld)))
DCM <- transform(DCM, id = as.numeric(factor(juld)))

#Modified DCM cleaning
for (i in 1:length(modified_dcm_juld_to_remove)){
  rho_DCM_modified_densityprofiledf <- rho_DCM_modified_densityprofiledf[!(rho_DCM_modified_densityprofiledf$juld == modified_dcm_juld_to_remove[i]),]
  modified_DCMprofiles <- modified_DCMprofiles[!(modified_DCMprofiles$juld == modified_dcm_juld_to_remove[i]),]
  modified_DCM <- modified_DCM[!(modified_DCM$juld == modified_dcm_juld_to_remove[i]),]
  }

rho_DCM_modified_densityprofiledf <- transform(rho_DCM_modified_densityprofiledf,id=as.numeric(factor(juld)))
modified_DCMprofiles <- transform(modified_DCMprofiles, id = as.numeric(factor(juld)))
modified_DCM <- transform(modified_DCM, id = as.numeric(factor(juld)))

# CHECK VISU
#Check Visu DCM
tmp <- DCMprofiles[DCMprofiles$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
off_fluo <- tmp$fluo[which.min(tmp$fluo[1:depthindex])]
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex]-off_fluo)
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, fgauss(x,tmp$Fsurf[1], tmp$Zdemi[1],
                        tmp$Fmax[1], tmp$Zmax[1], tmp$dz[1]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i + 1

#Check Visu modified DCM
tmp <- modified_DCMprofiles[modified_DCMprofiles$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
off_fluo <- tmp$fluo[which.min(tmp$fluo[1:depthindex])]
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$fluo[1:depthindex]-off_fluo)
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, fgauss(x,tmp$Fsurf[1], tmp$Zdemi[1],
                        tmp$Fmax[1], tmp$Zmax[1], tmp$dz[1]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i + 1


#CREATION OF TEXT FILES FOR FURTHER DIVA ANALYSIS
TOTAL_DCM_PROFILES <- rbind(DCMprofiles, modified_DCMprofiles)
TOTAL_DCM_PROFILES <- transform(TOTAL_DCM_PROFILES,id=as.numeric(factor(juld)))
TOTAL_DCM <- rbind(DCM, modified_DCM)
TOTAL_DCM <- transform(TOTAL_DCM, id = as.numeric(factor(juld)))
TOTAL_DCM <- TOTAL_DCM[order(TOTAL_DCM$id),]
TOTAL_RHO_PROFILES <- rbind(rho_DCM_densityprofiledf, rho_DCM_modified_densityprofiledf)
TOTAL_RHO_PROFILES <- transform(TOTAL_RHO_PROFILES, id = as.numeric(factor(juld)))
TOTAL_RHO_PROFILES <- TOTAL_RHO_PROFILES[order(TOTAL_RHO_PROFILES$id),]
TOTAL_RHO <- rbind(rhoDCM, rho_modified_DCM)
TOTAL_RHO <- transform(TOTAL_RHO, id = as.numeric(factor(juld)))
TOTAL_RHO <- TOTAL_RHO[order(TOTAL_RHO$id),]


TOTAL_DCM_PROFILES <- ddply(TOTAL_DCM_PROFILES,~month, transform, season=1*(month %in% c(12,1,2))+
                 2*(month %in% c(3,4,5 ))+
                 3*(month %in% c(6,7,8 ))+
                 4*(month %in% c(9,10,11 )))

TOTAL_DCM <- ddply(TOTAL_DCM,~month, transform, season=1*(month %in% c(12,1,2))+
                              2*(month %in% c(3,4,5 ))+
                              3*(month %in% c(6,7,8 ))+
                              4*(month %in% c(9,10,11 )))

TOTAL_RHO_PROFILES <- ddply(TOTAL_RHO_PROFILES,~month, transform, season=1*(month %in% c(12,1,2))+
                              2*(month %in% c(3,4,5 ))+
                              3*(month %in% c(6,7,8 ))+
                              4*(month %in% c(9,10,11 )))

TOTAL_RHO <- ddply(TOTAL_RHO,~month, transform, season=1*(month %in% c(12,1,2))+
                              2*(month %in% c(3,4,5 ))+
                              3*(month %in% c(6,7,8 ))+
                              4*(month %in% c(9,10,11 )))

#TEXT FILES creation FOR DIVA

# DCM ONLY

subDCM <- subset(DCM, select = c("lon", "lat", "juld", "Zmax","DOY", "year",
                                 "month", "day", "file"))
sub_rhoDCM <- subset(rhoDCM, select = c("lon", "lat", "juld", "rho_dcm","DOY", "year",
                                     "month", "day", "filename"))

dcm_rho_depth <- cbind(subDCM, sub_rhoDCM$rho_dcm)

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
sub_rho_modified_DCM <- subset(rho_modified_DCM, select = c("lon", "lat", "juld", "rho_dcm","DOY", "year",
                                        "month", "day", "filename"))

modified_dcm_rho_depth <- cbind(sub_modified_DCM, sub_rho_modified_DCM$rho_dcm)

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

sub_TOTAL_RHO <- subset(TOTAL_RHO, select = c("lon", "lat", "juld", "rho_dcm","DOY", "year",
                                                     "month", "day", "filename"))

TOTAL_DCM_DEPTH_RHO <- cbind(sub_TOTAL_DCM, sub_TOTAL_RHO$rho_dcm)

# write.table(sub_TOTAL_DCM, file="TOTAL_DCM_DEPTH.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(sub_TOTAL_RHO, file="TOTAL_DCM_RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# write.table(TOTAL_DCM_DEPTH_RHO, file="TOTAL_DCM_DEPTH&RHO.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)


############################################################
######################## GRAPHICS ##########################
############################################################

# Normalize data
normalize <- function(data){
  data <- (data - min(data, na.rm = T))/(max(data, na.rm = T) - min(data, na.rm = T))
}

# NORMALIZED DATASETS

#DCMprofiles

NORM_DCM <- normalize(DCMprofiles$fluo)
DCMprofiles2 <- cbind(DCMprofiles, NORM_DCM)

NORM_modified_DCM <- normalize(modified_DCMprofiles$fluo)
modified_DCMprofiles2 <- cbind(modified_DCMprofiles$fluo)

NORM_TOT_DCM <- normalize(TOTAL_DCM_PROFILES$fluo)
TOTAL_DCM_PROFILES2 <- cbind(TOTAL_DCM_PROFILES, NORM_TOT_DCM)

NORM_TOT_RHO <- normalize(TOTAL_RHO_PROFILES$FLUO_ADJUSTED)
TOTAL_RHO_PROFILES2 <- cbind(TOTAL_RHO_PROFILES, NORM_TOT_RHO)


# Skype Marilaure
# Faire les trucs de MLD min et Max (cfr Navarro) + regarder la PAR (1 et 10%)

#Density graphics -------


#Temporal evolution for 2017
tmp <- TOTAL_RHO_PROFILES2[TOTAL_RHO_PROFILES2$year == 2017,]
#tmp <- TOTAL_RHO_PROFILES2
Per1 <- tmp[tmp$DOY <= 30,]
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

#same as DepthDataExtraction
ggplot(perioddf, aes(x=density, y=NORM_TOT_RHO, group=juld)) +
  geom_point() + facet_grid(~group) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Fluorescence (rfu)") +
  coord_flip() + scale_x_reverse() #+ ylim(0,500)

#Seasonal Analysis per year

#Ex : 2017

season_2017_rho <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2017,]
ggplot(season_2017_rho, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  coord_flip() + scale_color_viridis() +  scale_x_reverse()

season_2016_rho <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2016,]
ggplot(season_2016_rho, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  coord_flip() + scale_color_viridis() +  scale_x_reverse() + ylim(0,400)

season_2015_rho <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2015,]
ggplot(season_2015_rho, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  coord_flip() + scale_color_viridis() +  scale_x_reverse() + ylim(0,400)

season_2014_rho <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2014,]
ggplot(season_2014_rho, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  coord_flip() + scale_color_viridis() +  scale_x_reverse() + ylim(0,400)

#Conclusion : Ne se focaliser que sur les années 2015 à 2017.

##### REMOVE MONTHS 10 to 1 ####

test <- densitygaussianprofilesdf
test <- test[!test$month == 10,]
test <- test[!test$month == 11,]
test <- test[!test$month == 12,]
test <- test[!test$month == 1,]
densitygaussianprofilesdf <- test
densitygaussianprofilesdf <- transform(densitygaussianprofilesdf,id=as.numeric(factor(juld)))

####

ggplot(densitygaussianprofilesdf, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_color_viridis() +  scale_x_reverse()

ggplot(densitygaussianprofilesdf, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_color_gradientn(colours = rainbow(5)) + scale_x_reverse()

#for continuous scale
ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=season, group=juld)) +
  geom_point() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_color_viridis() + scale_x_reverse()

#discrete case
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 1] <- "Winter"
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 2] <- "Spring"
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 3] <- "Summer"
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 4] <- "Autumn"

ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=season, group=juld)) +
  geom_point(size = 0.2) + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_x_reverse()

#NOTE : quid des mois hors max et min MLD (soit mois d'octobre à janvier)

densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 1] <- "Jan"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 2] <- "Feb"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 3] <- "Mar"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 4] <- "Apr"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 5] <- "May"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 6] <- "Jun"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 7] <- "Jul"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 8] <- "Aug"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 9] <- "Sep"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 10] <- "Oct"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 11] <- "Nov"
densitygaussianprofilesdf2$month[densitygaussianprofilesdf2$month == 12] <- "Dec"

densitygaussianprofilesdf2$month = factor(densitygaussianprofilesdf2$month, levels=c('Jan','Feb','Mar','Apr','May','Jun',
                                                                                     'Jul','Aug','Sep','Oct','Nov','Dec'))

ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_x_reverse()

#2017

tmp <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2017,]

ggplot(tmp, aes(x=density, y=chla, color=season, group=juld)) +
  geom_point(size = 0.1) + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_x_reverse()

ggplot(tmp, aes(x=density, y=chla, color=season, group=juld)) +
  geom_point() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_x_reverse()

ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,400) +
  coord_flip() + scale_color_gradientn(colours = rainbow(5)) + scale_x_reverse()
coord_flip() + scale_color_viridis() + scale_x_reverse()

#Depth graphics -------
#plots (cfr. Navarro et al. 2013, figure 3)

# ggplot(profiledf, aes(x=-depth, y=chla, color=month, group=juld)) + geom_line() +
#   facet_grid(~year) + xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-80, 0) + coord_flip() + scale_colour_gradientn(colours = terrain.colors(6))
# 
# ggplot(profiledf, aes(x=-depth, y=chla, color=month, group=juld)) + geom_line() +
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-80, 0) + coord_flip() + scale_colour_brewer()
# 
# ggplot(profiledf, aes(x=-depth, y=chla, color=as.factor(month), group=as.factor(juld))) +
#   geom_line() +
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-80, 0) + coord_flip() + scale_colour_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12))
# 
# ggplot(profiledf, aes(x=-depth, y=chla, color=as.factor(month), group=as.factor(juld))) +
#   geom_line() + facet_grid(~year)+
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-80, 0) + coord_flip() + scale_colour_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12))
# 
# ggplot(profiledf, aes(x=-depth, y=chla, color=as.factor(month), group=as.factor(juld))) +
#   geom_line() +
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-80, 0) + coord_flip() + scale_colour_manual(values=c(1,2,3,4,5,6,7,8,9,10,11,12))

#this one
ggplot(profiledf, aes(x=depth, y=chla, color=month, group=juld)) +
  geom_point() + facet_grid(~year) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) + ylim(0,400) +
  coord_flip() + scale_color_gradientn(colours = rainbow(5))

# ggplot(profiledf, aes(x=-depth, y=chla, color=month, group=juld)) +
#   geom_line() + facet_grid(~year) +
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-80, 0) + 
#   coord_flip() + scale_color_gradientn(colours = rainbow(5))

#or this one
ggplot(profiledf, aes(x=depth, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) +ylim(0,400) +
  coord_flip() + scale_color_viridis()

# ggplot(profiledf, aes(x=-depth, y=chla, color=month, group=juld)) +
#   geom_line() + 
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-50, -10) + coord_flip() + scale_color_gradientn(colours = rainbow(5))

#seasonal analysis (clearer view)
profiledf2<-ddply(profiledf,~month, transform, season=1*(month %in% c(12,1,2))+
                    2*(month %in% c(3,4,5 ))+
                    3*(month %in% c(6,7,8 ))+
                    4*(month %in% c(9,10,11 )))


# profiledf2$season[profiledf2$season == 1] <- "Winter"
# profiledf2$season[profiledf2$season == 2] <- "Spring"
# profiledf2$season[profiledf2$season == 3] <- "Summer"
# profiledf2$season[profiledf2$season == 4] <- "Autumn"


#OR THIS ONE ALSO (1 = winter)
ggplot(profiledf2, aes(x=depth, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) +ylim(0,400) +
  coord_flip() + scale_color_viridis()


#Temporal evolution for 2017
tmp <- profiledf2[profiledf2$year == 2017,]
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

ggplot(perioddf, aes(x=depth, y=chla, group=juld)) +
  geom_line() + facet_grid(~group) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) +ylim(0,400) +
  coord_flip() 

ggplot(perioddf, aes(x=depth, y=chla, group=juld)) +
  geom_point(size = 0.8) + facet_grid(~group) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) +ylim(0,400) +
  coord_flip() 

#Graphique Black Sea --------------
# Nasty code to get a user-defined season
argodf2<-ddply(argodf,~month, transform, season=1*(month %in% c(12,1,2))+
                 2*(month %in% c(3,4,5 ))+
                 3*(month %in% c(6,7,8 ))+
                 4*(month %in% c(9,10,11 )))

#pour tenir compte des saisons
# argodf4 <- subset(argodf2, select=c("lon","lat","juld","rho_anomalymax","season"))
# 

# write.table(argodf4, file="rhoInput_season.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

# argodf5 <- subset(argodf2, select=c("lon","lat","juld","rho_anomalymax","DOY"))
# write.table(argodf5, file="rhoInput_DOY.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

argodf2$season[which(argodf2$season == 1)]<-"Winter"
argodf2$season[which(argodf2$season == 2)]<-"Spring"
argodf2$season[which(argodf2$season == 3)]<-"Summer"
argodf2$season[which(argodf2$season == 4)]<-"Autumn"

# Le fait de le définir en "factor" va faire qu'il va respecter l'ordre imposé
argodf2$season<-factor(argodf2$season,levels = c("Winter","Spring","Summer","Autumn"))

# ggplot(argodf2, aes(x = rho_anomalymax, fill=factor(year))) + geom_density(alpha=.5) +
#   facet_grid(.~season,scales = "free") + coord_flip() + xlim(c(16,12))

#Spatial coverage of ARGO buoys ----------------------------------
#Map libraries
library(ggmap)
library(ggalt)

#Black Sea coordinates
bs <- c(26.5,40,43,46)
myMap <-get_map(location=bs, source="google", crop=FALSE)
ggmap(myMap) +
  geom_point(aes(x=lon, y=lat, color=factor(season) ),
             data = argodf2, alpha = .8)+facet_grid(year~Platform)

myMap <-get_map(location=bs, source="google", maptype = "satellite", crop=FALSE)
ggmap(myMap) +
  geom_point(aes(x=lon, y=lat, color=factor(Platform)),
             data = argodf2, alpha = .8)


#DIVAND files creation ------------------

diva_DCM_depth <- subset(gaussdata, select=c("lon","lat","juld","Zmax"))
#Remove NA's
diva_DCM_depth <- diva_DCM_depth[complete.cases(diva_DCM_depth),]
write.table(diva_DCM_depth, file="diva_DCM_depth.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

diva_DCM_rho <- subset(rho_dcm, select=c("lon","lat","juld","rho_dcm"))
#Remove NA's
diva_DCM_rho <- diva_DCM_rho[complete.cases(diva_DCM_rho),]
write.table(diva_DCM_rho, file="diva_DCM_rho.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

#2014
tmp <- rho_dcm[rho_dcm$year == 2014,]
diva_DCM_rho_2014 <- subset(tmp, select=c("lon","lat","juld","rho_dcm"))
#Remove NA's
diva_DCM_rho_2014 <- diva_DCM_rho_2014[complete.cases(diva_DCM_rho_2014),]
write.table(diva_DCM_rho_2014, file="diva_DCM_rho_2014.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

#2015
tmp <- rho_dcm[rho_dcm$year == 2015,]
diva_DCM_rho_2015 <- subset(tmp, select=c("lon","lat","juld","rho_dcm"))
#Remove NA's
diva_DCM_rho_2015 <- diva_DCM_rho_2015[complete.cases(diva_DCM_rho_2015),]
write.table(diva_DCM_rho_2015, file="diva_DCM_rho_2015.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

#2016
tmp <- rho_dcm[rho_dcm$year == 2016,]
diva_DCM_rho_2016 <- subset(tmp, select=c("lon","lat","juld","rho_dcm"))
#Remove NA's
diva_DCM_rho_2016 <- diva_DCM_rho_2016[complete.cases(diva_DCM_rho_2016),]
write.table(diva_DCM_rho_2016, file="diva_DCM_rho_2016.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

#2017
tmp <- rho_dcm[rho_dcm$year == 2017,]
diva_DCM_rho_2017 <- subset(tmp, select=c("lon","lat","juld","rho_dcm"))
#Remove NA's
diva_DCM_rho_2017 <- diva_DCM_rho_2017[complete.cases(diva_DCM_rho_2017),]
write.table(diva_DCM_rho_2017, file="diva_DCM_rho_2017.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

#Some statistics  --------------

#Monthly repartition before gaussian profiles were retained
ggplot(gaussiandf, aes(month)) + geom_histogram(binwidth = 0.5) + 
  labs(title="All profiles [908]")

#Monthly repartition after non gaussian profiles were eliminated
ggplot(gaussdata, aes(month)) + geom_histogram(binwidth = 0.5) + 
  labs(title="Only gaussian profiles [668]")

#Introduction of seasons
# "Nasty-code" from Arthur
# Nasty code to get a user-defined season
gaussdata_season <- ddply(gaussdata,~month, transform, season=1*(month %in% c(12,1,2))+
                            2*(month %in% c(3,4,5 ))+
                            3*(month %in% c(6,7,8 ))+
                            4*(month %in% c(9,10,11 )))

rho_dcm_season <- ddply(rho_dcm,~month, transform, season=1*(month %in% c(12,1,2))+
                          2*(month %in% c(3,4,5 ))+
                          3*(month %in% c(6,7,8 ))+
                          4*(month %in% c(9,10,11 )))

#Seasonal variability for depth and density, for each season compared to all data (all years)
#Depth
mean_depth_winter <- mean(gaussdata_season[gaussdata_season$season == 1,]$Zmax, na.rm = T)
sd_depth_winter <- sd(gaussdata_season[gaussdata_season$season == 1,]$Zmax, na.rm = T)
mean_depth_spring <- mean(gaussdata_season[gaussdata_season$season == 2,]$Zmax, na.rm = T)
sd_depth_spring <- sd(gaussdata_season[gaussdata_season$season == 2,]$Zmax, na.rm = T)
mean_depth_summer <- mean(gaussdata_season[gaussdata_season$season == 3,]$Zmax, na.rm = T)
sd_depth_summer <- sd(gaussdata_season[gaussdata_season$season == 3,]$Zmax, na.rm = T)
mean_depth_autumn <- mean(gaussdata_season[gaussdata_season$season == 4,]$Zmax, na.rm = T)
sd_depth_autumn <- sd(gaussdata_season[gaussdata_season$season == 4,]$Zmax, na.rm = T)

#Density anomaly
mean_rho_winter <- mean(rho_dcm_season[rho_dcm_season$season == 1,]$rho_dcm, na.rm = T)
sd_rho_winter <- sd(rho_dcm_season[rho_dcm_season$season == 1,]$rho_dcm, na.rm = T)
mean_rho_spring <- mean(rho_dcm_season[rho_dcm_season$season == 2,]$rho_dcm, na.rm = T)
sd_rho_spring <- sd(rho_dcm_season[rho_dcm_season$season == 2,]$rho_dcm, na.rm = T)
mean_rho_summer <- mean(rho_dcm_season[rho_dcm_season$season == 3,]$rho_dcm, na.rm = T)
sd_rho_summer <- sd(rho_dcm_season[rho_dcm_season$season == 3,]$rho_dcm, na.rm = T)
mean_rho_autumn <- mean(rho_dcm_season[rho_dcm_season$season == 4,]$rho_dcm, na.rm = T)
sd_rho_autumn <- sd(rho_dcm_season[rho_dcm_season$season == 4,]$rho_dcm, na.rm = T)

#All years
mean_depth_all <- mean(gaussdata$Zmax, na.rm =T)
sd_depth_all <- sd(gaussdata$Zmax, na.rm =T)
mean_rho_all <- mean(rho_dcm$rho_dcm, na.rm = T)
sd_rho_all <- sd(rho_dcm$rho_dcm, na.rm = T)

#Comparison between rho_max_MLD (from Luc) VS rho_max_ARGO [DCM] 
#YEAR 2015
mean_rho_2015 <- mean(rho_dcm[rho_dcm$year == 2015,]$rho_dcm, na.rm =T)
sd_rho_2015 <- sd(rho_dcm[rho_dcm$year == 2015,]$rho_dcm, na.rm =T)

#YEAR 2016
mean_rho_2016 <- mean(rho_dcm[rho_dcm$year == 2016,]$rho_dcm, na.rm =T)
sd_rho_2016 <- sd(rho_dcm[rho_dcm$year == 2016,]$rho_dcm, na.rm =T)

mean_rho_2014 <- mean(rho_dcm[rho_dcm$year == 2014,]$rho_dcm, na.rm =T)
sd_rho_2014 <- sd(rho_dcm[rho_dcm$year == 2014,]$rho_dcm, na.rm =T)

mean_rho_2017 <- mean(rho_dcm[rho_dcm$year == 2017,]$rho_dcm, na.rm =T)
sd_rho_2017 <- sd(rho_dcm[rho_dcm$year == 2017,]$rho_dcm, na.rm =T)


########################################################################################
########################################################################################

# DATA INTERCALIBRATION ################################################################

########################################################################################
########################################################################################

# Only 3 floats -> WOD 6900807, 6901866 and the now one 6902340 (will soon be on the DAC/Coriolis)
# Comparison firstly based on temperature and salinity

# Get salinity and temperature profiles from all BCG-floats equipped with CDOM sensor

#Get each profile from ARGO data -----
GET_PROFILES <- function(file, VAR_NAME){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
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
  
  ### Direct use of adjusted values if available
  vardf <- ExtractVar(paste0(VAR_NAME,"_ADJUSTED"),FloatInfo)
  if (all(is.na(vardf)) == TRUE) {
    vardf <- ExtractVar(VAR_NAME,FloatInfo)
  }
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  #Construction du data frame final a lieu ici
  data.frame(depth           = vardf$depth,
             juld            = vardf$juld,
             var             = vardf$value,
             qc              = vardf$qc,
             day             = month.day.year(vardf$juld,c(1,1,1950))$day,
             month           = month.day.year(vardf$juld,c(1,1,1950))$month,
             year            = month.day.year(vardf$juld,c(1,1,1950))$year,
             DOY             = as.integer(strftime(as.Date(vardf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon             = vardf$lon,
             lat             = vardf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
}

TEMP_profiles <- ldply(as.list(filename), function(file){
  temp <- GET_PROFILES(file,"TEMP")
})
colnames(TEMP_profiles)[3] <- "temp"
colnames(TEMP_profiles)[4] <- "qc_temp"

PSAL_profiles <- ldply(as.list(filename), function(file){
  psal <- GET_PROFILES(file,"PSAL")
})
colnames(PSAL_profiles)[3] <- "psal"
colnames(PSAL_profiles)[4] <- "qc_psal"

TS_profiles <- TEMP_profiles
TS_profiles$psal <- PSAL_profiles$psal
TS_profiles$qc_psal <- PSAL_profiles$qc_psal

TS_profiles <- transform(TS_profiles,id=as.numeric(factor(juld)))

ggplot(TS_profiles, aes(x=depth, y=temp, group=juld)) +
  geom_line() + 
  xlab("Depth (m)") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(TS_profiles, aes(x=depth, y=psal, group=juld)) +
  geom_line() + 
  xlab("Depth (m)") + ylab("Salinity") +
  coord_flip() + scale_x_reverse()

#remove data with QC = 3 and QC = 4
#TO DO

CHLA_profiles <- ldply(as.list(filename), function(file){
  temp <- GET_PROFILES(file,"CHLA")
})
colnames(CHLA_profiles)[3] <- "chla"

ggplot(CHLA_profiles, aes(x=depth, y=chla, group=juld)) +
  geom_line() + 
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  coord_flip() + scale_x_reverse() + geom_hline(yintercept=0)

# Database of profiles co-located in space and in time
# LAUCH COORDINATES
deployment_lon <- 28.88 #29 CTD
deployment_lat <- 43.09 #43.10 CTD

# carré de 0.1° de côté ~ 10 km de côté
# Database of profiles co-located in space 

#PAR RAPPORT À L'ENDROIT DU DÉPLOIEMENT
co_space_deployment <- TS_profiles[TS_profiles$lat <= 43.2 & TS_profiles$lat >= 43
                                   & TS_profiles$lon <= 29.2 & TS_profiles$lon >= 29,]

#remove rows containing NA's
co_space_deployment <- co_space_deployment[complete.cases(co_space_deployment),]

#remove duplicates
co_space_deployment<- unique(co_space_deployment)

ggplot(co_space_deployment, aes(x=depth, y=temp, color = id, group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(co_space_deployment, aes(x=depth, y=temp, color = factor(id), group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(co_space_deployment, aes(x=depth, y=psal, color = id, group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Salinity") +
  coord_flip() + scale_x_reverse()

ggplot(co_space_deployment, aes(x=depth, y=psal, color = factor(id), group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Salinity") +
  coord_flip() + scale_x_reverse()

#TS diagram
ggplot(co_space_deployment, aes(x=psal, y=temp, color = factor(id), group=juld)) +
  geom_point() + 
  xlab("Salinity") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()


# co_space_deployment <- TS_profiles[TS_profiles$lat <= 43.3 & TS_profiles$lat >= 42.9
#                                    & TS_profiles$lon <= 29.4 & TS_profiles$lon >= 28.8,]


# Database of profiles co-located in time

co_time_deployment <- TS_profiles[TS_profiles$DOY <= 93 & 
                                    TS_profiles$DOY >= 83,]

#TS diagram
ggplot(co_time_deployment, aes(x=psal, y=temp, color = id, group=juld)) +
  geom_point() + 
  xlab("Salinity") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(co_time_deployment, aes(x=psal, y=temp, color = factor(id), group=juld)) +
  geom_point() + 
  xlab("Salinity") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(co_time_deployment, aes(x=depth, y=temp, color = id, group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(co_time_deployment, aes(x=depth, y=temp, color = factor(id), group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Temperature (°C)") +
  coord_flip() + scale_x_reverse()

ggplot(co_time_deployment, aes(x=depth, y=psal, color = id, group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Salinity") +
  coord_flip() + scale_x_reverse()

ggplot(co_time_deployment, aes(x=depth, y=psal, color = factor(id), group=juld)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Salinity") +
  coord_flip() + scale_x_reverse()

###############
path = "/home/flo/R/TFE/tfe/6903240"
new_float <- list.files(path = path, pattern="*.nc")

#plus proche du déploiement disponible jusqu'ici...

filetest <- new_float[1]

#Opening the file in a open-only mode
ncfile   <<- nc_open(filetest, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)

juld     <- ncvar_get(ncfile,"JULD")
pres     <- as.data.frame(ncvar_get(ncfile,"PRES"))
lon      <- ncvar_get(ncfile,"LONGITUDE")
lat      <- ncvar_get(ncfile,"LATITUDE")
chla_adjusted     <- as.data.frame(ncvar_get(ncfile,"CHLA_ADJUSTED"))
chla     <- as.data.frame(ncvar_get(ncfile,"CHLA")) 
psal_adjusted     <- as.data.frame(ncvar_get(ncfile,"PSAL_ADJUSTED"))
psal     <- as.data.frame(ncvar_get(ncfile,"PSAL"))
temp_adjusted     <- as.data.frame(ncvar_get(ncfile,"TEMP_ADJUSTED"))
temp     <- as.data.frame(ncvar_get(ncfile,"TEMP"))
cdom_adjusted <- as.data.frame(ncvar_get(ncfile,"CDOM_ADJUSTED"))
cdom <- as.data.frame(ncvar_get(ncfile,"CDOM"))

bsdf <- cbind(pres, chla_adjusted, chla, psal_adjusted, psal, 
              temp_adjusted, temp, cdom_adjusted, cdom, juld, lon, lat)

colnames(bsdf) <- c("depth","chla_adjusted","chla","psal_adjusted",
                    "psal","temp_adjusted","temp","cdom_adjusted",
                    "cdom", "juld", "lon", "lat")

ggplot(bsdf, aes(x=depth, y=cdom_adjusted)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Temperature") +
  coord_flip() + scale_x_reverse()

#En attendant que le fichier soit sur le DAC tout beau tout clean
#Essai de "clean" de données

temp_index <- which(bsdf$temp > 10)
psal_index <- which(bsdf$psal > 25)

bsdf$temp[temp_index] <- NA
bsdf$psal[psal_index] <- NA

#repérage indice chla grâce au cdom (à l'air logique quand on 
#analyse les chiffres)
cdom_index <- which(bsdf$cdom_adjusted < 1)

bsdf$cdom_adjusted[cdom_index] <- NA
bsdf$chla[cdom_index] <- NA

a <- ggplot(bsdf, aes(x=depth, y=chla)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Chloro") +
  coord_flip() + scale_x_reverse()

#NPQ correction for GOLD profile

psalGOLD <- gsw_SA_from_SP(bsdf$psal,bsdf$depth,bsdf$lon,bsdf$lat)
tempGOLD <- gsw_CT_from_t(psalGOLD,bsdf$temp,bsdf$depth)
rho_anomalyGOLD <- gsw_sigma0(psalGOLD,tempGOLD)
t7 <- as.data.frame(cbind(rho_anomalyGOLD, bsdf$depth))

#MLD est de +/- 20-21 mètres donc ça colle avec le "faux" DCM que l'on voit (données argo
#pas encore traitées)

t8 <-which(bsdf$depth < 20 & bsdf$depth > 0)

b <- ggplot(bsdf[t8,], aes(x=depth, y=chla)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Chloro") +
  coord_flip() + scale_x_reverse()


grid.arrange(a,b,ncol=2, nrow = 1)

#################
# Function who plots T and S from profiles chosen (either manually
# or automatically)
#

ggplot(bsdf, aes(x=depth)) + 
  geom_point(aes(y = temp, colour = "Temperature")) +
  geom_point(aes(y = psal, colour = "Salinity")) + coord_flip() +
  xlab("Depth (m)") + ylab("Temperature (°C)") + 
  scale_x_reverse()

bad_data <- which(co_space_deployment$qc_temp == 4 
                  | co_space_deployment$qc_temp == 3)  
co_space_deployment$temp[bad_data] <- NA

trace <- subset(co_space_deployment[co_space_deployment$id == 173,], select = c("depth","temp"))

ggplot(trace, aes(x=depth, y=temp)) +
  geom_point() + 
  xlab("Depth (m)") + ylab("Temperature") +
  coord_flip() + scale_x_reverse() + geom_point(data = trace)

############# SEARCH OF GOOD PROFILES IN 100KM AND APRIL-MARS ####

#space
search_launch <- TS_profiles[TS_profiles$lat <= 44 & TS_profiles$lat >= 42
                             & TS_profiles$lon <= 30.1 & TS_profiles$lon >= 28.1,]

#time (AVRIL-MARS PAS DISPO...) --> février et mai.. 
search_launch <- search_launch[search_launch$month == 2 |
                                 search_launch$month == 5,]

#REMOVE QC 3 AND 4
bad_data <- which(search_launch$qc_temp == 4 
                  | search_launch$qc_temp ==3)
search_launch$temp[bad_data] <- NA

bad_data <- which(search_launch$qc_psal == 4 
                  | search_launch$qc_psal ==3)
search_launch$psal[bad_data] <- NA

search_launch <- transform(search_launch,id=as.numeric(factor(juld)))

i <- 6

tmp <- search_launch[search_launch$id == i,]
tmp <- subset(tmp, select = c("depth","temp", "psal","juld"))
bsdf2 <- subset(bsdf, select = c("depth","temp","psal","juld"))
compare <- rbind(tmp,bsdf2)

a <- ggplot(compare, aes(x=depth, y=temp, colour = factor(juld))) +
  geom_point(size = 0.5) + 
  xlab("Depth (m)") + ylab("Temperature") +
  coord_flip() + scale_x_reverse() 

b <- ggplot(compare, aes(x=depth, y=psal, colour =factor(juld))) +
  geom_point(size = 0.5) + 
  xlab("Depth (m)") + ylab("Salinity") +
  coord_flip() + scale_x_reverse() 

chla_test <- CHLA_profiles[CHLA_profiles$juld == tmp$juld[1],]
chla_test <- subset(chla_test, select = c("depth","chla", "juld"))
bsdf3 <- subset(bsdf, select = c("depth", "chla", "juld"))

test <- rbind(chla_test, bsdf3)

c <- ggplot(test, aes(x=depth, y=chla, colour = factor(juld))) +
  geom_point(size = 0.5) + 
  xlab("Depth (m)") + ylab("Chla") +
  coord_flip() + scale_x_reverse() 

grid.arrange(a,b,c,ncol=2, nrow = 2)

i <- i + 1