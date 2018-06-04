##############################################################################################################
#                                                   REMARKS                                                  #
# THE STRUCTURE OF EXTRACT_VAR HAD TO BE ADAPTED BECAUSE WE DO NOT DEAL WITH MERGED ARGO FILES BUT INSTEAD   #
# WE WORK WITH INDIVIDUAL PROFILE FILES. MOREOVER, THE STRUCTURE OF INDIVIDUAL FILES WAS NOT UNIFORM         #
#                                                                                                            #
# AIM OF THE SCRIPT : INVESTIGATION OF MLD PROPERTIES IN THE BLACK SEA SINCE 2005 (ARGO DATA BEGAN)          #
##############################################################################################################

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
library(pbapply)

#get working directory
wd <- getwd()
filename <- list.files(path = paste0(wd,"/Flo_floats"), pattern="*.nc")

#get variable
ExtractVar <- function(Var,FloatInfo){
  with(FloatInfo,{
    lvar             <- as.vector(ncvar_get(ncfile,Var))
    lvar_qc          <- as.vector(ncvar_get(ncfile,paste0(Var,"_QC")))
    lvar_qctab       <- llply(lvar_qc,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab<-do.call(cbind,lvar_qctab)
    lvar_qctab <- as.vector(lvar_qctab)
    lvar_qc_pres     <- as.vector(ncvar_get(ncfile,paste0("PRES","_QC")))
    lvar_qctab_pres  <- llply(lvar_qc_pres,function(qcstring){
      as.numeric(unlist(strsplit(qcstring,split="")))
    })
    lvar_qctab_pres <- do.call(cbind,lvar_qctab_pres)
    lvar_qctab_pres <- as.vector(lvar_qctab_pres)
   
    data.frame(value    = lvar,
               qc       = lvar_qctab,
               qc_pres       = lvar_qctab_pres,
               depth    = pres,
               variable = Var,
               juld     = juld,
               lon      = lon,
               lat      = lat)
    })

  }

#get density profiles
density_profiles <- ldply(as.list(filename),function(file){
  
  ncfile   <<- nc_open(paste0(wd,"/Flo_floats/",file), write = FALSE, verbose = FALSE, suppress_dimvals = FALSE)
  juld     <- mean(ncvar_get(ncfile,"JULD"))
  pres     <- as.vector(ncvar_get(ncfile,"PRES"))
  lon      <- mean(ncvar_get(ncfile,"LONGITUDE"))
  lat      <- mean(ncvar_get(ncfile,"LATITUDE"))
  FloatInfo<-list(juld = mean(juld),
                  pres=pres,
                  lon=mean(lon),
                  lat=mean(lat))
  
  tempdf <- ExtractVar("TEMP",FloatInfo)
  tempdf$value[which(tempdf$qc == 4 | tempdf$qc == 3)] <- NA
  psaldf <- ExtractVar("PSAL",FloatInfo)
  psaldf$value[which(psaldf$qc == 4 | psaldf$qc == 3)] <- NA
  
  nc_close(ncfile)
  
  psal <- gsw_SA_from_SP(psaldf$value,psaldf$depth,psaldf$lon,psaldf$lat)
  temp <- gsw_CT_from_t(psal,tempdf$value,tempdf$depth)
  sigma <- gsw_sigma0(psal,temp)

  show(file)
  
  data.frame(sigma   = sigma,
             depth   = tempdf$depth, 
             juld    = tempdf$juld, 
             qc_pres = tempdf$qc_pres,
             lon     = tempdf$lon,
             lat     = tempdf$lat,
             day     = month.day.year(tempdf$juld,c(1,1,1950))$day,
             month   = month.day.year(tempdf$juld,c(1,1,1950))$month,
             year    = month.day.year(tempdf$juld,c(1,1,1950))$year,
             DOY     = as.integer(strftime(as.Date(tempdf$juld,origin = '1950-01-01'), 
                                         format ="%j")))#Day Of Year
  
  
})

#REMOVE WRONG PROFILES & clean
save <- density_profiles
density_profiles <- density_profiles[!(density_profiles$lon < 27.5),]
density_profiles <- density_profiles[!rowSums(is.na(density_profiles[,1:2]))==2,]
density_profiles <- density_profiles[!(density_profiles$qc_pres == 4),]
rownames(density_profiles) <- NULL

#density profiles summary
argodf<- ddply(density_profiles,~juld,summarize,
                    check = length(which(depth <= 10)),
                    bottomdepth = max(depth),
                    min_depth = min(depth),
                    lon=mean(lon),
                    lat=mean(lat),
                    day = day[1], month = month[1],
                    year = year[1], DOY = DOY[1])

crit1 <- argodf[(argodf$min_depth > 10),]#all are stuck profiles except one
crit2 <- argodf[(argodf$bottomdepth < 300),]
crit3 <- argodf[(argodf$check <2),]
criteria <- unique(rbind(crit1, crit2, crit3))

#remove bad profiles
for (i in 1:length(criteria$juld)){
  argodf <- argodf[!(argodf$juld == criteria$juld[i]),]
  density_profiles <- density_profiles[!(density_profiles$juld == criteria$juld[i]),]
}

#reorder
density_profiles <- transform(density_profiles,id=as.numeric(factor(juld)))
rownames(argodf) <- NULL
rownames(density_profiles) <- NULL

#MLD crteria
sigma_criteria <- 0.03
depth_ref <- 10

#MLD computation
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
  
  show(i)
  
  data.frame(sigma_surface = sigma_surface, MLD = MLD,
             juld = tmp$juld[1], lon = tmp$lon[1], lat = tmp$lat[1], day = month.day.year(tmp$juld[1],c(1,1,1950))$day,
             month = month.day.year(tmp$juld[1],c(1,1,1950))$month,
             year = month.day.year(tmp$juld[1],c(1,1,1950))$year,
             DOY = as.integer(strftime(as.Date(tmp$juld[1],origin = '1950-01-01'), 
                                       format ="%j")))
})

#reorder
MLDdf <- MLDdf[!rowSums(is.na(MLDdf[,1:2]))==2,]
rownames(MLDdf) <- NULL

#get readable time
new_time2 <- as.Date(MLDdf$juld, origin = as.Date("1950-01-01"))
MLDdf$date <- new_time2

#Black Sea coordinates
bs <- c(26.5,40,43,46)
myMap <-get_map(location=bs, source="google", crop=FALSE, maptype = "terrain")
ggmap(myMap) +
  geom_point(aes(x=lon, y=lat, color=factor(month)), size = 1,
             data = MLDdf, alpha = 1) + #+facet_grid(file) + scale_color_manual(values=c("red", "green", "yellow","blue"))
labs(x = "Longitude (°)", y = "Latitude (°)", color = "Month\n") 

years <- unique(MLDdf$year)

# Some stats  
MLD_sigma <- ddply(tmp,~month~year,summarize,
            maxMLD = max(MLD),
            sigma_surface = sigma_surface[which.max(MLD)],
            juld = juld[which.max(MLD)])

#readable time
new_time <- as.Date(MLD_sigma$juld, origin = as.Date("1950-01-01"))
MLD_sigma$date <- new_time

ggplot(MLD_sigma, aes(x = new_time, y = maxMLD, color = factor(month))) + 
  geom_point() + geom_hline(yintercept = 14) 

MLDdf <- MLDdf[!(MLDdf$MLD > 100),]

ggplot(MLDdf, aes(x = new_time2, y = MLD, color = factor(month))) + 
  geom_point() + geom_hline(yintercept = 14) 

meanMLD_month <- ddply(MLDdf, ~month, summarize,
                       meanMLD = mean(MLD, na.rm = T),
                       sdMLD= sd(MLD, na.rm = T),
                       data = length(MLD),
                       maxMLD = max(MLD, na.rm = T),
                       sigmasurf = mean(sigma_surface, na.rm = T))
