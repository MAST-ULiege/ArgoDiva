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

#CORRECTION OF XING ET AL. 2017
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


cdom_floats <- c("6900807_Mprof.nc","6901866_Mprof.nc")


cdom_profiles <- ldply(as.list(cdom_floats),function(file){
  
  #Opening the file in a open-only mode
  ncfile   <<- nc_open(file, write = FALSE, verbose = FALSE, suppress_dimvals = FALSE)
  
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
  
  
  cdomdf <- ExtractVar("CDOM", FloatInfo) 
  chladf <- ExtractVar("CHLA", FloatInfo)
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  show(file)
  
  data.frame(depth    = cdomdf$depth,
             juld     = cdomdf$juld,
             fluo     = chladf$value,
             cdom     = cdomdf$value,
             qc_cdom  = cdomdf$qc,
             qc_chla  = chladf$qc,
             day      = month.day.year(cdomdf$juld,c(1,1,1950))$day,
             month    = month.day.year(cdomdf$juld,c(1,1,1950))$month,
             year     = month.day.year(cdomdf$juld,c(1,1,1950))$year,
             DOY      = as.integer(strftime(as.Date(cdomdf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon      = cdomdf$lon,
             lat      = cdomdf$lat,
             dir      = cdomdf$dir,
             Platform = as.numeric(unique(id)),
             type     = "Argo")
})

cdom_profiles <- cdom_profiles[-(which(cdom_profiles$qc_chla == 4)),] 
cdom_profiles <- cdom_profiles[-(which(cdom_profiles$dir == "D")),]

#ROW NAMES != ROW NUMBERS !!!! ---> REORDER !!
#cdom_profiles <- cdom_profiles[order(cdom_profiles$juld),]
rownames(cdom_profiles) <- NULL

mmed <- function(x,n=5){runmed(x,n)}

smoothed_fluo <- ldply(as.list(unique(cdom_profiles$juld)), function(i){
  tmp <- cdom_profiles[cdom_profiles$juld == i,]
  tmp <- mmed(tmp$fluo, 5)
  data.frame(smoothed_fluo = tmp)
})

smoothed_cdom <- ldply(as.list(unique(cdom_profiles$juld)), function(i){
  tmp <- cdom_profiles[cdom_profiles$juld == i,]
  tmp <- mmed(tmp$cdom, 5)
  data.frame(smoothed_cdom = tmp)
})

#REPLACE FLUO WITH SMOOTHED FLUO & CDOM
cdom_profiles[,3] <- smoothed_fluo
cdom_profiles[,4] <- smoothed_cdom

cdom_argodf<- ddply(cdom_profiles,~juld,summarize,
                    qc_chla = qc_chla[which.max(fluo)],
                    depthmax = depth[which.max(fluo)],
                    maxvalue = fluo[which.max(fluo)],
                    depthmin = depth[which.max(fluo):length(fluo)][which.min(fluo[which.max(fluo):length(fluo)])],
                    integration = sum(fluo),
                    bottomdepth = max(depth),
                    min_depth = min(depth),
                    dir = dir[1],
                    lon=mean(lon),
                    lat=mean(lat),
                    Platform = Platform[1],
                    day = day[1], month = month[1],
                    year = year[1], DOY = DOY[1])

cdom_criteriadf1 <- cdom_argodf[(cdom_argodf$min_depth > 5),]#all stuck profiles except one
cdom_criteriadf2 <- cdom_argodf[(cdom_argodf$bottomdepth < 800),]
cdom_criteriadf <- unique(rbind(cdom_criteriadf1, cdom_criteriadf2))

#remove profiles that have been categorized as bad profiles
for (i in 1:length(cdom_criteriadf$juld)){
  cdom_argodf <- cdom_argodf[!(cdom_argodf$juld == cdom_criteriadf$juld[i]),]
  cdom_profiles <- cdom_profiles[!(cdom_profiles$juld == cdom_criteriadf$juld[i]),]
}

rownames(cdom_profiles) <- NULL
rownames(cdom_argodf) <- NULL
cdom_profiles <- transform(cdom_profiles,id=as.numeric(factor(juld)))
cdom_argodf <- transform(cdom_argodf,id=as.numeric(factor(juld)))

#FOR CHECK MLD
cdom_density_profiles <- ldply(as.list(cdom_floats),function(file){
  
  ncfile   <<- nc_open(file, write = FALSE, verbose = FALSE, suppress_dimvals = FALSE)
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
  
  #DENSITY ANOMALY BASED ON TEOS-10
  psal <- gsw_SA_from_SP(psaldf$value,psaldf$depth,psaldf$lon,psaldf$lat)
  temp <- gsw_CT_from_t(psal,tempdf$value,tempdf$depth)
  sigma <- gsw_sigma0(psal,temp)
  
  show(file)
  
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
                                         format ="%j")))#Day Of Year
  
})

good_profiles <- cdom_argodf$juld

cdom_density_profiles <- cdom_density_profiles[cdom_density_profiles$juld %in% good_profiles,]
rownames(cdom_density_profiles) <- NULL
cdom_density_profiles <- transform(cdom_density_profiles,id=as.numeric(factor(juld)))

sigma_criteria <- 0.03
depth_ref <- 10


MLDdf <- ldply(as.list(1:length(unique(cdom_density_profiles$id))), function(i){
  
  tmp <- cdom_density_profiles[cdom_density_profiles$id == i,]
  
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

#EACH PROFILES HAS BEEN VALIDATED SO FAR ==> CALIBRATION EXERCISE STARTS NOW

#FDOM-based method (Xing et al. 2017)

FDOMdf <- ldply(as.list(1:length(unique(cdom_density_profiles$id))), function(i){
  
  tmp <- cdom_profiles[cdom_profiles$id == i,]
  MaxDepth <- cdom_argodf$bottomdepth[i]#Max depth of the profile
  TopDepth <- cdom_argodf$depthmin[i]#Apparent minimum 
  
  #Compute correlation between FDOM and FChla (to test the hypothesis
  #that FDOM perturbates linearly FChla at depth)
  
  #Linear regression between TopDepth and MaxDepth
  #FCHLA_MEASURED = SLOPE_FDOM * FDOM_MEASURED + C
  
  calibrange <- tmp[which(tmp$depth == TopDepth):which(tmp$depth == MaxDepth),]
  linearMod <- lm(fluo ~ cdom, data = calibrange)
  slope_fdom <- coef(linearMod)[[2]]
  C <- coef(linearMod)[[1]]
  #summary(linearMod)
  
  data.frame(slope_fdom = slope_fdom, C = C)
})

# chla_estimation <- slope_fdom*calibrange$cdom + C
# 
# r = 1-sum((calibrange$fluo - chla_estimation)^2)/(sum((calibrange$fluo - mean(calibrange$fluo))^2))
# 
# scatterplot(fluo ~ cdom, data = calibrange)
# 
# ggplot(calibration, aes(x=depth, y=chla, colour = origin)) +
#   geom_point(size = 1) + 
#   xlab("Depth (m)") + ylab("CHLA") +
#   coord_flip() + scale_x_reverse() 


#DESCENTE PROFILE DURING DEPLOYMENT
maxparam <- max(count.fields("ogsbio008a_001_01_T250.csv", sep = '')) 
launch_param <- read.csv(file = "ogsbio008a_001_01_T250.csv", sep = '', 
                         na.strings="NA", 
                         col.names = paste0("V", seq_len(maxparam), 
                                            fill=TRUE), header = FALSE, stringsAsFactors = FALSE)

maxparam <- max(count.fields("ogsbio008a_001_00_05.csv", sep = '')) 
launch_descent <- read.csv(file = "ogsbio008a_001_00_05.csv", sep = '', 
                           na.strings="NA", 
                           col.names = paste0("V", seq_len(maxparam), 
                                              fill=TRUE), header = TRUE,
                           stringsAsFactors = FALSE)

launch_descent <- subset(launch_descent,
                         select = c("V1TRUE","V2TRUE","V8TRUE","V9TRUE",
                                    "V16TRUE","V17TRUE","V19TRUE"))

colnames(launch_descent) <- c("depth","date","temp","psal","PAR","chla","cdom")

launch_descent <- launch_descent[which(launch_descent$date > 2.018e+13),]
launch_descent <- launch_descent[!rowSums(is.na(launch_descent[,2:7]))==5,]
launch_descent <- launch_descent[launch_descent$depth >= 0,]

dark_chla <- as.numeric(launch_param$V38TRUE[5])
scale_chla <- as.numeric(launch_param$V37TRUE[5])
launch_descent$chla <- (launch_descent$chla - dark_chla)*scale_chla
dark_cdom <- as.numeric(launch_param$V42TRUE[5])
scale_cdom <- as.numeric(launch_param$V41TRUE[5])
launch_descent$cdom <- (launch_descent$cdom - dark_cdom)*scale_cdom

descent <- subset(launch_descent, select = c("depth","temp","psal","chla","cdom"))
descent <- descent[!rowSums(is.na(descent[,4:5]))==2,]


maxparam <- max(count.fields("ogsbio008a_001_00_09.csv", sep = '')) 
launch_ascent <- read.csv(file = "ogsbio008a_001_00_09.csv", sep = '', 
                          na.strings="NA", 
                          col.names = paste0("V", seq_len(maxparam), 
                                             fill=TRUE), header = TRUE,
                          stringsAsFactors = FALSE)

launch_ascent <- subset(launch_ascent,
                        select = c("V1TRUE","V2TRUE","V8TRUE","V9TRUE",
                                   "V16TRUE","V17TRUE","V19TRUE"))

colnames(launch_ascent) <- c("depth","date","temp","psal","PAR","chla","cdom")

launch_ascent <- launch_ascent[which(launch_ascent$date > 2.018e+13),]
launch_ascent <- launch_ascent[!rowSums(is.na(launch_ascent[,2:7]))==5,]
launch_ascent <- launch_ascent[launch_ascent$depth >= 0,]

dark_chla <- as.numeric(launch_param$V38TRUE[5])
scale_chla <- as.numeric(launch_param$V37TRUE[5])
launch_ascent$chla <- (launch_ascent$chla - dark_chla)*scale_chla
dark_cdom <- as.numeric(launch_param$V42TRUE[5])
scale_cdom <- as.numeric(launch_param$V41TRUE[5])
launch_ascent$cdom <- (launch_ascent$cdom - dark_cdom)*scale_cdom

ascent <- subset(launch_ascent, select = c("depth","temp","psal","chla","cdom"))
ascent$origin <- "ASCENT"
ascent <- subset(ascent, select = c("depth","temp","psal","chla","cdom"))
ascent <- ascent[!rowSums(is.na(ascent[,2:3]))==2,]
# write.table(descent, file="descent_launch_cleaned.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

rownames(descent) <- NULL
rownames(ascent) <- NULL



ggplot(descent, aes(x = chla, y = depth)) + geom_path() +
  scale_y_reverse()

ggplot(descent, aes(x = cdom, y = depth)) + geom_path() +
  scale_y_reverse()

#Spike test (QC test 9) 
#initialize empty array (test <- ()) si on modifie le code un jour
chla <- descent$chla
spike <- sapply(3:(length(chla)-3),function(i){
  Test_Value <- abs(chla[i]-((chla[i+1]+chla[i-1])/2)) - abs((chla[i+1]-chla[i-1])/2)
  Threshold_Value <- abs(median(c(chla[i-2], chla[i-1], chla[i], chla[i+1], chla[i+2]), na.rm=TRUE)) +
    abs(sd(c(chla[i-2], chla[i-1], chla[i], chla[i+1], chla[i+2]), na.rm=TRUE))
  if (Test_Value > Threshold_Value){
    index1 <- i
  }
})

index1 <- unlist(spike)

#Spike test (version 2014 (+ 2016 pour meilleure explication)) PAS LE MEME !! NEGATIVE spike test
#positive spikes are considered as OK (non flaggée mais je peux vous assurer qu'il le faudrait sur
#certaines)
RES <- sapply(3:(length(chla)-3),function(i){
  #plus de valeurs absolues?
  tmp <- chla[i] - median(c(chla[i-2], chla[i-1], chla[i], chla[i+1], chla[i+2]), na.rm=TRUE)
})

#RES <- sort(as.vector(RES), decreasing = F)
Q10 <- quantile(RES, probs = 0.1)

spike2 <- sapply(1:(length(RES)),function(i){
  if(RES[i] < 2*Q10)
    index4 <- i + 2
})  

index4 <- unlist(spike2)

index <- c(index1, index4)

#remove potentially bad value
descent <- descent[-index,]

var_average <- function(df){
  
  tmp <- df[order(df$depth),]
  tmp <- tmp[!rowSums(is.na(tmp[,4:5])) == 2,]
  tmp2 <- aggregate(chla ~ depth, data = tmp, FUN = mean)
  tmp2b <- aggregate(chla ~ depth, data = tmp, FUN = sd, na.rm = T)
  tmp3 <- aggregate(cdom ~ depth, data = tmp, FUN = mean)
  tmp3b <- aggregate(cdom ~ depth, data = tmp, FUN = sd, na.rm = T)
  duplicate <- which(duplicated(tmp$depth))
  data.frame(depth = tmp2$depth, chla = tmp2$chla, 
             chla_sd = tmp2b$chla, cdom = tmp3$cdom, 
             cdom_sd = tmp3b$cdom)
}

var_average_density <- function(df){
  
  tmp <- df[order(df$depth),]
  tmp <- tmp[!rowSums(is.na(tmp[,2:3])) == 2,]
  tmp2 <- aggregate(temp ~ depth, data = tmp, FUN = mean)
  tmp2b <- aggregate(temp ~ depth, data = tmp, FUN = sd, na.rm = T)
  tmp3 <- aggregate(psal ~ depth, data = tmp, FUN = mean)
  tmp3b <- aggregate(psal ~ depth, data = tmp, FUN = sd, na.rm = T)
  duplicate <- which(duplicated(tmp$depth))
  data.frame(depth = tmp2$depth, temp = tmp2$temp, 
             temp_sd = tmp2b$temp, psal = tmp3$psal, 
             psal_sd = tmp3b$psal)
}

descent <- var_average(descent)
ascent <- var_average_density(ascent)

psal <- gsw_SA_from_SP(ascent$psal,ascent$depth,29,43.1)
temp <- gsw_CT_from_t(psal,ascent$temp, ascent$depth)
sigma <- gsw_sigma0(psal,temp)
sigma_surface <- approx(ascent$depth,sigma,10)$y
ascent <- cbind(ascent, sigma)
MLD <- max(ascent$depth[ascent$sigma <= (sigma_surface + sigma_criteria)], na.rm = T)

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

descent$cdom <- mmed(descent$cdom, 5)
descent$chla <- mmed(descent$chla, 5)

ggplot(descent, aes(x = chla, y = depth)) + geom_path() +
  scale_y_reverse()

ggplot(descent, aes(x = cdom, y = depth)) + geom_path() +
  scale_y_reverse()

#250m
rownames(descent) <- NULL
MaxDepth <- max(descent$depth)#Max depth of the profile
TopDepth <- descent$depth[which.min(descent$chla)]

#Compute correlation between FDOM and FChla (to test the hypothesis
#that FDOM perturbates linearly FChla at depth)

#Linear regression between TopDepth and MaxDepth
#FCHLA_MEASURED = SLOPE_FDOM * FDOM_MEASURED + C

calibrange <- descent[which(descent$depth == TopDepth):which(descent$depth == MaxDepth),]
linearMod <- lm(chla ~ cdom, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]

chla_cor1 <- descent$chla - (slope_fdom*descent$cdom) - C
descent <- cbind(descent, chla_cor1)
chla_cor_npq <- quenching_correction(descent$chla_cor1, descent$depth, MLD)
chla_cor2 <- chla_cor_npq
descent <- cbind(descent, chla_cor2)
chla_cor3 <- descent$chla - (mean(FDOMdf$slope_fdom) * descent$cdom) - mean(FDOMdf$C)
descent <- cbind(descent, chla_cor3)

depth_hplc <- c(0,30,50,70,90,100,140,200,250)
chla_hplc <- c(0.6111,0.6433,0.3090,0.0675,0.0229,0.0168,0.012,0.0032,
               0.0035)
hplc <- as.data.frame(cbind(depth_hplc, chla_hplc))

ggplot(descent, aes(x = chla, y = depth)) + geom_path(size = 0.5) +
  scale_y_reverse() + geom_path(data = descent, aes(x = chla_cor1, y = depth), colour = "blue", size = 0.5) +
  geom_point(data=hplc, aes(y = depth_hplc, x = chla_hplc), colour = "red", size = 2, shape = 15) +
  xlab(expression(Chlorophyll~a~(mg/m^3)))+ ylab("Depth (m)") + geom_path(data = descent,
                                                                          aes(x = chla_cor3, y = depth),
                                                                          colour = "green", size = 0.5) +
  theme(text=element_text(size=12))

indexdepth <- vector()

for (i in 1:length(hplc$chla)){
  indexdepth[i] <- which.min(abs(descent$depth - hplc$depth_hplc[i]))
}

rmse_chla_vec <- descent$chla[indexdepth]
rmse_chla_vec <- descent$chla_cor1[indexdepth]
rmse_chla_vec <- descent$chla_cor2[indexdepth]
rmse_chla_vec <- descent$chla_cor3[indexdepth]


rmse <- sqrt(sum((hplc$chla_hplc - rmse_chla_vec)^2)/length(hplc$depth_hplc))

#RMSE for top layer correction
rmse <- sqrt(sum((hplc$chla_hplc[1:6] - rmse_chla_vec[1:6])^2)/length(hplc$depth_hplc[1:6]))

#for deep
rmse <- sqrt(sum((hplc$chla_hplc[7:9] - rmse_chla_vec[7:9])^2)/length(hplc$depth_hplc[7:9]))

#below 50m
rmse <- sqrt(sum((hplc$chla_hplc[4:9] - rmse_chla_vec[4:9])^2)/length(hplc$depth_hplc[4:9]))

# #NPQ correction
# 
# 
# 
# 
# #HPLC
# # depth_hplc <- c(0,30,50,70,90,100,140,200,250,500,750,1000)
# # chla_hplc <- c(0.6111,0.6433,0.3090,0.0675,0.0229,0.0168,0.012,0.0032,
# #                0.0035,0.0023,0.0039,0.0042)
# depth_hplc <- c(0,30,50,70,90,100,140,200,250)
# chla_hplc <- c(0.6111,0.6433,0.3090,0.0675,0.0229,0.0168,0.012,0.0032,
#                0.0035)
# temp_hplc <- rep(NA, length(depth_hplc))
# psal_hplc <- rep(NA, length(depth_hplc))
# cdom_hplc <- rep(NA, length(depth_hplc))
# 
# hplcprofile <- as.data.frame(cbind(depth_hplc, temp_hplc, psal_hplc,
#                                    chla_hplc, cdom_hplc))
# hplcprofile$origin <- "HPLC"
# colnames(hplcprofile) <- c("depth","temp","psal","chla","cdom","origin")
# 
# ggplot(descent_avg, aes(x = chla, y = depth)) +
#   geom_point() + scale_y_reverse() +
#   geom_point(data = chla_cor, aes(x = chla, y = depth, color = "red"))+
#   geom_point(data = hplcprofile, aes(x = chla, y = depth, color = "yellow"))
# 
# indexdepth <- vector()
# 
# for (i in 1:length(hplcprofile$chla)){
# indexdepth[i] <- which.min(abs(chla_cor$depth - hplcprofile$depth[i]))
# }
# 
# rmse_chla_vec <- chla_cor$chla[indexdepth]
# rmse_chla_vec <- descent_avg$chla[indexdepth]

# rmse = sqrt(sum((hplcprofile$chla - rmse_chla_vec)^2)/length(hplcprofile$depth))
# 
# 
# #ASCENT PROFILE 1000 M
# path = "/home/flo/R/TFE/tfe/6903240"
# new_float <- list.files(path = path, pattern="*.nc")
# 
# #plus proche du déploiement disponible jusqu'ici...
# 
# filetest <- new_float[7]
# 
# #Opening the file in a open-only mode
# ncfile   <<- nc_open(filetest, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
# 
# juld     <- ncvar_get(ncfile,"JULD")
# pres     <- as.data.frame(ncvar_get(ncfile,"PRES"))
# lon      <- ncvar_get(ncfile,"LONGITUDE")
# lat      <- ncvar_get(ncfile,"LATITUDE")
# chla     <- as.data.frame(ncvar_get(ncfile,"CHLA"))
# psal     <- as.data.frame(ncvar_get(ncfile,"PSAL"))
# temp     <- as.data.frame(ncvar_get(ncfile,"TEMP"))
# cdom <- as.data.frame(ncvar_get(ncfile,"CDOM"))
# 
# GOLD <- cbind(pres, chla, psal, temp, cdom, juld, lon, lat)
# 
# colnames(GOLD) <- c("depth","chla",
#                     "psal","temp",
#                     "cdom", "juld", "lon", "lat")
# 
# GOLD$cdom[GOLD$cdom == 99999] <- NA
# GOLD <- GOLD[!rowSums(is.na(GOLD[,2:5])) == 4,]
# GOLD <- GOLD[!(is.na(GOLD[,3]) == FALSE),]
# 
# #APPLICATION OF CORRECTION BASED ON THE MEAN OF 404 PROFILES FROM THE BLACK SEA
# chla_cor <- as.data.frame(GOLD$chla - (0.035*GOLD$cdom) - mean(FDOMdf$C))/2
# colnames(chla_cor) <- "chla"
# chla_cor$depth <- GOLD$depth
# 
# #HPLC
# depth_hplc <- c(0,30,50,70,90,100,140,200,250,500,750,1000)
# chla_hplc <- c(0.6111,0.6433,0.3090,0.0675,0.0229,0.0168,0.012,0.0032,
#                0.0035,0.0023,0.0039,0.0042)
# temp_hplc <- rep(NA, length(depth_hplc))
# psal_hplc <- rep(NA, length(depth_hplc))
# cdom_hplc <- rep(NA, length(depth_hplc))
# 
# hplcprofile <- as.data.frame(cbind(depth_hplc, temp_hplc, psal_hplc, 
#                                    chla_hplc, cdom_hplc))
# hplcprofile$origin <- "HPLC"
# colnames(hplcprofile) <- c("depth","temp","psal","chla","cdom","origin")
# 
# ggplot(GOLD, aes(x = chla, y = depth)) +
#   geom_point() + scale_y_reverse() + 
#   geom_point(data = chla_cor, aes(x = chla, y = depth, color = "red"))+
#   geom_point(data = hplcprofile, aes(x = chla, y = depth, color = "yellow"))
# 
# for (i in 1:length(hplcprofile$chla)){
#   indexdepth[i] <- which.min(abs(chla_cor$depth - hplcprofile$depth[i]))
# }
# 
# rmse_chla_vec <- chla_cor$chla[indexdepth]
# rmse_chla_vec <- GOLD$chla[indexdepth]
# 
# hplcprofile[1:3,] <- NA
# rmse_chla_vec[1:3] <- NA
# 
# rmse <- sqrt(sum((hplcprofile$chla - rmse_chla_vec)^2, na.rm = T)/
#               (length(hplcprofile$depth) - length(which(is.na(hplcprofile$chla) == TRUE))))
# 


#CDOM analysis
# tmp <- cdom_profiles[cdom_profiles$id == i,]
# ggplot(tmp, aes(x = cdom, y = depth)) + geom_point() + 
#   scale_y_reverse()
# i <- i + 1


cdom_profiles<- ddply(cdom_profiles,~month, transform, season=1*(month %in% c(12,1,2))+
                        2*(month %in% c(3,4,5 )) +
                        3*(month %in% c(6,7,8 )) +
                        4*(month %in% c(9,10,11 )))

a <- ggplot(cdom_profiles, aes(x = cdom, y =depth, colour = factor(Platform))) + 
  geom_point(size = 0.1) + scale_y_reverse() + 
  xlab(expression(CDOM~(ppb)))+ ylab("Depth (m)") + theme(legend.position = "none") +
  theme(text=element_text(size=12)) + scale_color_manual(values = c("steelblue4", "steelblue2"))
  #scale_x_continuous(breaks = c(0, 0.5))
  # theme(axis.title.x=element_blank(), axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())


b <- ggplot(cdom_profiles, aes(x = fluo, y =depth, colour = factor(Platform))) + 
  geom_point(size = 0.1) + scale_y_reverse() + xlim(c(0,2)) +
  xlab(expression(Chlorophyll~a~(mg/m^3)))+ ylab("Depth (m)") + theme(legend.position = "none") +
  theme(text=element_text(size=12)) + scale_color_manual(values = c("green4", "green2")) + 
#scale_x_continuous(breaks = c(0, 0.5))
theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
      axis.ticks.y=element_blank())

grid.arrange(a,b, ncol=2, nrow = 1)

#Creation of a basin scale CDOM profile (no real seasonality has been spotted, at least in the deeper part)
bindf <- ldply(as.list(1:length(unique(cdom_profiles$id))), function(i){
  tmp <- cdom_profiles[cdom_profiles$id == i,]
  tmp1 <- tmp[tmp$depth <= 200,]
  tmp2 <- tmp[tmp$depth > 200 & tmp$depth <= 400,]
  #tmp3 <- tmp[tmp$depth > 400 & tmp$depth <= 1000,]#on peut allonger un peu loin (regarder bottom max des profils sans CDOM)
  tmp3 <- tmp[tmp$depth > 400,]
  bin1 <- as.data.frame(binAverage(tmp1$depth, tmp1$cdom, xmin = 0, xmax = floor(max(tmp1$depth, na.rm = T)), xinc = 1))
  bin2 <- as.data.frame(binAverage(tmp2$depth, tmp2$cdom, xmin = 200, xmax = floor(max(tmp2$depth, na.rm = T)), xinc = 5))
  bin3 <- as.data.frame(binAverage(tmp3$depth, tmp3$cdom, xmin = 400, xmax = floor(max(tmp3$depth, na.rm = T)), xinc = 10))
  bin <- rbind(bin1, bin2, bin3)
  data.frame(depth = bin$x, cdom = bin$y, juld = tmp$juld[1], platform = tmp$Platform[1])
})

bin2df <- ddply(bindf, ~depth, summarize,
                #l = length(which(!is.na(depth) == TRUE)),
                cdom = mean(cdom, na.rm = T))

ggplot(bin2df, aes(x = cdom, y =depth)) + 
  geom_point(size = 0.5) + scale_y_reverse() + geom_smooth(span=1)

# ggplot(bin2df, aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=depth<<-..y..))
# 
# ggplot(bin2df, aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=cdom<<-..x..))
# 
# ggplot(cdom_profiles[cdom_profiles$Platform == 6901866,], aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=cdom<<-..x..))
# 
# ggplot(cdom_profiles[cdom_profiles$Platform == 6900807,], aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=cdom1<<-..x..))
# 
# ggplot(cdom_profiles[cdom_profiles$Platform == 6900807,], aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=depth1<<-..y..))
# 
# ggplot(cdom_profiles[cdom_profiles$Platform == 6901866,], aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=cdom2<<-..x..))
# 
# ggplot(cdom_profiles[cdom_profiles$Platform == 6901866,], aes(x = cdom, y =depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=depth2<<-..y..))
# 
# ggplot(bin2df, aes(x = cdom, y = depth)) + 
#   geom_point(size = 0.5) + scale_y_reverse() + geom_smooth(span=1)
# 
# cdom <- bin2df$cdom
# depth <- bin2df$depth
# plot(x=cdom, y = -depth)
# fit3 <- lm(-depth~poly(cdom, 3, raw=T))
# fit4 <- lm(-depth ~ exp(cdom))
# xx <- seq(1, 15, length = 50)
# lines(xx, predict(fit3, data.frame(cdom = xx)), col="red")
# lines(xx, predict(fit4, data.frame(cdom = xx)), col="blue")

#CHOICE TO BE MADE -> bin ONLY one platform CDOM profiles (for consistency at the bottom)
cdom_test <- cdom_profiles[cdom_profiles$Platform == 6900807,]
cdom_test <- transform(cdom_test,id=as.numeric(factor(juld)))
c <- ggplot(cdom_test, aes(x = cdom, y =depth)) + 
  geom_point(size = 0.5) + scale_y_reverse() + stat_smooth(span=1, aes(outfit=depth<<-..y..)) +
  ylab("Depth (m)") + xlab("CDOM (ppb)") + theme(text=element_text(size=12))  

bindf3 <- ldply(as.list(1:length(unique(cdom_test$id))), function(i){
  tmp <- cdom_test[cdom_test$id == i,]
  tmp1 <- tmp[tmp$depth <= 200,]
  tmp2 <- tmp[tmp$depth > 200 & tmp$depth <= 400,]
  #tmp3 <- tmp[tmp$depth > 400 & tmp$depth <= 1000,]#on peut allonger un peu loin (regarder bottom max des profils sans CDOM)
  tmp3 <- tmp[tmp$depth > 400,]
  bin1 <- as.data.frame(binAverage(tmp1$depth, tmp1$cdom, xmin = 0, xmax = floor(max(tmp1$depth, na.rm = T)), xinc = 1))
  bin2 <- as.data.frame(binAverage(tmp2$depth, tmp2$cdom, xmin = 200, xmax = floor(max(tmp2$depth, na.rm = T)), xinc = 5))
  bin3 <- as.data.frame(binAverage(tmp3$depth, tmp3$cdom, xmin = 400, xmax = floor(max(tmp3$depth, na.rm = T)), xinc = 10))
  bin <- rbind(bin1, bin2, bin3)
  data.frame(depth = bin$x, cdom = bin$y, juld = tmp$juld[1], platform = tmp$Platform[1])
})

bin3df <- ddply(bindf3, ~depth, summarize,
                #l = length(which(!is.na(depth) == TRUE)),
                cdom = mean(cdom, na.rm = T))

d <- ggplot(bin3df, aes(x = cdom, y = depth)) + 
  geom_point(size = 1) + scale_y_reverse() + #geom_smooth(span=1) +
  xlab("CDOM (ppb)") + theme(text=element_text(size=12)) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

grid.arrange(c,d,nrow = 1, ncol = 2)

#better
ggplot(cdom_test, aes(x = cdom, y =depth)) + 
  geom_point(size = 0.5) + scale_y_reverse() + #stat_smooth(span=1, aes(outfit=depth<<-..y..)) +
  ylab("Depth (m)") + xlab("CDOM (ppb)") + theme(text=element_text(size=12)) +
  geom_point(data = bin3df, aes(x = cdom, y = depth), colour = "red", size = 0.5)


# cdom <- bin3df$cdom
# depth <- bin3df$depth
# plot(x=cdom, y = -depth)
# fit3 <- lm(-depth~poly(cdom, 5, raw=T))
# fit4 <- lm(-depth ~ exp(cdom))
# xx <- seq(1, 15, length = 50)
# lines(xx, predict(fit3, data.frame(cdom = xx)), col="red")
# lines(xx, predict(fit4, data.frame(cdom = xx)), col="blue")
                                                           
#1ST TEST : NO EQUATION FIT, only get bin values (true only mean profile)
cdom_validation <- cdom_profiles[cdom_profiles$Platform == 6901866,]
cdom_validation <- transform(cdom_validation,id=as.numeric(factor(juld)))

tmp <- cdom_validation[cdom_validation$id == i,]
tmp$fluo <- mmed(tmp$fluo, 5)
ggplot(tmp, aes(x = fluo, y =depth)) + 
  geom_path(size = 0.5) + scale_y_reverse() 

MaxDepth <- max(tmp$depth)#Max depth of the profile
#TopDepth <- tmp$depth[which.min(tmp$fluo)]
TopDepth <- tmp$depth[which.max(tmp$fluo):length(tmp$fluo)][which.min(tmp$fluo[which.max(tmp$fluo):length(tmp$fluo)])]

test <- data.frame(depth = approx(bin3df$depth, bin3df$cdom, tmp$depth)$x, 
                   cdom = approx(bin3df$depth, bin3df$cdom, tmp$depth)$y)

tmp$cdom_virtual <- test$cdom

calibrange <- tmp[which(tmp$depth == TopDepth):which(tmp$depth == MaxDepth),]

linearMod <- lm(fluo ~ cdom_virtual, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]

fluo_test_cor <- tmp$fluo - (slope_fdom*tmp$cdom_virtual) - C
tmp$fluo_virtual <- fluo_test_cor

ggplot(tmp, aes(x = fluo, y = depth)) + geom_point() + scale_y_reverse() + 
  geom_point(data = tmp, aes(x = fluo_virtual, y = depth), colour = "red") +
  geom_vline(xintercept = 0, lty = "dashed")

#i <- i + 1

#VISU : on voit que ça marche plutôt pas mal
#TO DO : calculer la correction sur le vrai profil de CDOM sur le flotteur 6901866
#FAIRE LA FORMULE DE RMSE ET ENSUITE RECALCULER UNE RMSE TOTALE, BELOW DE CHLA MINIMUM, ABOVE THE CHLA MINIUM
#THEN A-ONCE THIS METHO HAS BEEN CONSIDERED AS ACCURATE
# CORRECT NON CDOM FLOATS PROFILES (SHOW SOME GRAPHES)

#calibrange <- tmp[which(tmp$depth == TopDepth):length(tmp$depth),]
linearMod <- lm(fluo ~ cdom, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]
fluo_cor <- tmp$fluo - slope_fdom*tmp$cdom - C
tmp$true_fluo_cor <- fluo_cor

ggplot(tmp, aes(x = fluo, y = depth)) + geom_point() + scale_y_reverse() +
  geom_point(data = tmp, aes(x = fluo_virtual, y = depth), colour = "red") +
  geom_vline(xintercept = 0, lty = "dashed") + geom_point(data = tmp, aes(x = true_fluo_cor), colour = "green") + 
  xlab(expression(Chlorophyll~a~(mg/m^3))) + ylab("Depth (m)") +
    theme(text=element_text(size=12)) 


rmse_index <- which(!is.na(tmp$cdom_virtual == TRUE))
rmse <- sqrt(sum((tmp$true_fluo_cor[rmse_index] - tmp$fluo_virtual[rmse_index])^2)/length(tmp$fluo_virtual[rmse_index]))

#i == 88 

i <- i +  1

#data frame resulting from this analysis

cdom_virtual_validation <- ldply(as.list(1:length(unique(cdom_validation$id))), function(i){
  tmp <- cdom_validation[cdom_validation$id == i,]
  tmp$fluo <- mmed(tmp$fluo, 5)#note que je pourrais tout simplement faire cela avant pour TOUT le cdom_validation
  
  MaxDepth <- max(tmp$depth)#Max depth of the profile
  TopDepth <- tmp$depth[which.max(tmp$fluo):length(tmp$fluo)][which.min(tmp$fluo[which.max(tmp$fluo):length(tmp$fluo)])]
  
  #bin3df must be good defined
  test <- data.frame(depth = approx(bin3df$depth, bin3df$cdom, tmp$depth)$x, 
                     cdom = approx(bin3df$depth, bin3df$cdom, tmp$depth)$y)
  
  tmp$cdom_virtual <- test$cdom
  
  calibrange <- tmp[which(tmp$depth == TopDepth):which(tmp$depth == MaxDepth),]
  
  linearMod_virtual <- lm(fluo ~ cdom_virtual, data = calibrange)
  slope_fdom_virtual <- coef(linearMod_virtual)[[2]]
  C_virtual <- coef(linearMod_virtual)[[1]]
  
  fluo_test_cor <- tmp$fluo - (slope_fdom_virtual*tmp$cdom_virtual) - C_virtual
  tmp$fluo_virtual <- fluo_test_cor
  
  linearMod <- lm(fluo ~ cdom, data = calibrange)
  slope_fdom <- coef(linearMod)[[2]]
  C <- coef(linearMod)[[1]]
  fluo_cor <- tmp$fluo - slope_fdom*tmp$cdom - C
  tmp$true_fluo_cor <- fluo_cor
  
  ggplot(tmp, aes(x = fluo, y = depth)) + geom_point() + scale_y_reverse() +
    geom_point(data = tmp, aes(x = fluo_virtual, y = depth), colour = "red") +
    geom_vline(xintercept = 0, lty = "dashed") + geom_point(data = tmp, aes(x = true_fluo_cor), colour = "green")

  rmse_index <- which(!is.na(tmp$cdom_virtual == TRUE))
  rmse <- sqrt(sum((tmp$true_fluo_cor[rmse_index] - tmp$fluo_virtual[rmse_index])^2)/length(tmp$fluo_virtual[rmse_index]))
  
  data.frame(juld = tmp$juld[1], id = tmp$id[1], slope_true = slope_fdom, slope_virtual = slope_fdom_virtual,
             C_true = C, C_virtual = C_virtual, rmse = rmse)
  
})

#TEST OF VALIDATION WITH HPLC
#NOTE : on peut définir un profil moyen de plein de façon différentes
#en soi on pourrait également tenter d'inverser le profil de correction (ici 6900807) et de validation (6901866)

MaxDepth <- max(descent$depth)#Max depth of the profile
TopDepth <- descent$depth[which.min(descent$chla)]

#Compute correlation between FDOM and FChla (to test the hypothesis
#that FDOM perturbates linearly FChla at depth)

#Linear regression between TopDepth and MaxDepth
#FCHLA_MEASURED = SLOPE_FDOM * FDOM_MEASURED + C

calibrange <- descent[which(descent$depth == TopDepth):which(descent$depth == MaxDepth),]
linearMod <- lm(chla ~ cdom, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]


test <- data.frame(depth = approx(bin3df$depth, bin3df$cdom, descent$depth)$x, 
                   cdom = approx(bin3df$depth, bin3df$cdom, descent$depth)$y)

descent$cdom_virtual <- test$cdom

linearMod_virtual <- lm(chla ~ cdom_virtual, data = calibrange)
slope_fdom_virtual <- coef(linearMod_virtual)[[2]]
C_virtual <- coef(linearMod_virtual)[[1]]

fluo_test_cor <- descent$chla - (slope_fdom_virtual*descent$cdom_virtual) - C_virtual
descent$fluo_virtual <- fluo_test_cor

rmse <- sqrt(sum((descent$chla_cor1 - descent$fluo_virtual)^2)/length(tmp$fluo_virtual))

ggplot(descent, aes(x = chla, y = depth)) + geom_path(size = 1) +
  scale_y_reverse() + geom_path(data = descent, aes(x = chla_cor1, y = depth), colour = "blue", size = 1) +
  geom_point(data=hplc, aes(y = depth_hplc, x = chla_hplc), colour = "red", size = 2, shape = 15) +
  xlab(expression(Chlorophyll~a~(mg/m^3)))+ ylab("Depth (m)") + geom_path(data = descent,
                                                                          aes(x = fluo_virtual, y = depth),
                                                                          colour = "green", size = 1) +
  theme(text=element_text(size=12)) + ylim(c(100,0))

#cfr déf de indexdepth above

for (i in 1:length(hplc$chla)){
  indexdepth[i] <- which.min(abs(descent$depth - hplc$depth_hplc[i]))
}


rmse_chla_virtual_vec <- descent$fluo_virtual[indexdepth]

rmse <- sqrt(sum((hplc$chla_hplc - rmse_chla_virtual_vec)^2)/length(hplc$depth_hplc))


########## ALL PROFILES #################

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

##################

#test
nocdom <- profiles[profiles$Platform == 7900591 | profiles$Platform == 7900592,]
#nocdom <- profiles[profiles$Platform == 6900807 | profiles$Platform == 6901866,]
nocdom <- transform(nocdom,id=as.numeric(factor(juld)))
rownames(nocdom) <- NULL

#SMOOTHING FOR THEIR PROFILES !!
tmp <- nocdom[nocdom$id == i,]
tmp$fluo <- mmed(tmp$fluo, 5)
ggplot(tmp, aes(x = fluo, y =depth)) + 
  geom_path(size = 0.5) + scale_y_reverse() 

MaxDepth <- max(tmp$depth)#Max depth of the profile
#TopDepth <- tmp$depth[which.min(tmp$fluo)]
TopDepth <- tmp$depth[which.max(tmp$fluo):length(tmp$fluo)][which.min(tmp$fluo[which.max(tmp$fluo):length(tmp$fluo)])]

# test <- data.frame(depth = approx(bin2df$depth, bin2df$cdom, tmp$depth)$x,
#                    cdom = approx(bin2df$depth, bin2df$cdom, tmp$depth)$y)

test <- data.frame(depth = approx(bin3df$depth, bin3df$cdom, tmp$depth)$x,
                   cdom = approx(bin3df$depth, bin3df$cdom, tmp$depth)$y)

tmp$cdom <- test$cdom

calibrange <- tmp[which(tmp$depth == TopDepth):which(tmp$depth == MaxDepth),]

linearMod <- lm(fluo ~ cdom, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]

fluo_test_cor <- tmp$fluo - (slope_fdom*tmp$cdom) - C
tmp$fluo_cor <- fluo_test_cor

ggplot(tmp, aes(x = fluo, y = depth)) + geom_point() + scale_y_reverse() + 
  geom_point(data = tmp, aes(x = fluo_test_cor, y = depth), colour = "red") +
  geom_vline(xintercept = 0, lty = "dashed") + 
  xlab(expression(Chlorophyll~a~(mg/m^3)))+ ylab("Depth (m)")  +
  theme(text=element_text(size=16)) 
  #scale_x_continuous(breaks = c(0, 0.5))
  # theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank()) +
  # theme(axis.title.y=element_blank(), axis.text.y=element_blank(),
  #       axis.ticks.y=element_blank())


i <- i + 1


# #RESULTS -> profiles 108, 130,165, 167, 201, 212 are not good, everything else is OK
# check_bad_correction <- nocdom[nocdom$id %in% c(108,130,165, 167, 201, 212),]
# check_bad_correction <- ddply(check_bad_correction, ~juld, summarize,
#                               day = day[1],
#                               year = year[1],
#                               month = month[1],
#                               lat = lat[1],
#                               lon = lon[1],
#                               platform = Platform[1],
#                               surface = min(depth),
#                               bottom = max(depth),
#                               fluo_surface = fluo[1],
#                               depth_null = depth[which.min(fluo)],
#                               fluo_null = fluo[which.min(fluo)])

#Issue from the min depth (needs to be taken below the fluo MAX) -> note qu'en corrigeant d'abord le quenching, cette
# précaution n'est plus nécessaire mais le quenching correction doit se faire après Xing correction je dirais.. 

#####NOTE ADDITIONNELLE ###############################################################
#The mean profile stops at 950m so it is normal to non have any points after correction.
#######################################################################################
#

#######################################################
#SENSIBILITY TEST OF THE METHOD OF THE MEAN PROFILE####
#######################################################

#Test of the mean profile on the GOLD profile (even if CDOM is provided)
colnames(descent)[2] <- "fluo"
ggplot(descent, aes(x = fluo, y =depth)) + 
  geom_path(size = 0.5) + scale_y_reverse()
MaxDepth <- max(descent$depth)#Max depth of the profile
TopDepth <- descent$depth[which.max(descent$fluo):length(descent$fluo)][which.min(descent$fluo[which.max(descent$fluo):length(descent$fluo)])]

test <- data.frame(depth = approx(bin2df$depth, bin2df$cdom, descent$depth)$x, 
                   cdom = approx(bin2df$depth, bin2df$cdom, descent$depth)$y)

descent2 <- descent

descent2$virtualcdom <- test$cdom

calibrange <- descent2[which(descent2$depth == TopDepth):which(descent2$depth == MaxDepth),]

linearMod <- lm(fluo ~ virtualcdom, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]

fluo_cor <- descent2$fluo - (slope_fdom*descent2$virtualcdom) - C
descent2$fluo_cor <- fluo_cor

ggplot(descent2, aes(x = fluo, y = depth)) + geom_point() + scale_y_reverse() + 
  geom_point(data = descent2, aes(x = fluo_cor, y = depth), colour = "red") +
  geom_vline(xintercept = 0, lty = "dashed")

ggplot(descent2, aes(x = fluo, y = depth)) + geom_path() + scale_y_reverse() + 
  geom_path(data = descent2, aes(x = fluo_cor, y = depth), colour = "blue") +
  geom_vline(xintercept = 0, lty = "dashed") + 
  geom_point(data=hplc, aes(y = depth_hplc, x = chla_hplc),
             colour = "red", size = 2, shape = 15) + 
  geom_path(data = descent2, aes(x = chla_cor1, y = depth), colour = "green") 


######

#TRUE SENSITIVITY TEST#

#% Error =  |(T-E)/T|*100 where T = Corrected GOLD profile (with TRUE cdom profile)
#provided by the float

obs <- data.frame(depth = descent2$depth, fluo = descent2$chla_cor1, cdom = descent2$cdom)

#TEST 1 : use of mean cdom profile
exp <- bin2df

#extrapolation mean profile sur gold profile (only extrapolation on depth)
exp <- data.frame(depth = approx(exp$depth, exp$cdom, obs$depth)$x, 
                  cdom = approx(exp$depth, exp$cdom, obs$depth)$y)

#error_in = difference in % between real CDOM and virtual CDOM (mean CDOM profile)
#error_in = data.frame(depth = obs$depth, error = abs((obs$cdom-exp$cdom)/obs$cdom)*100)

rownames(descent2) <- NULL
MaxDepth <- max(descent2$depth)#Max depth of the profile
TopDepth <- descent2$depth[which.min(descent2$fluo)]

exp$fluo <- descent2$fluo

#Noise 
pc <- 0.5# %of noise added to the cdom signal (mean profile of cdom)
corrupt <- rbinom(length(exp$cdom),1,1)    # choose an average of 10% to corrupt at random
corrupt <- as.logical(corrupt)
noise <- rnorm(sum(corrupt),exp$cdom[corrupt],pc) # generate the noise to add
exp$noisy_cdom <- noise 
# exp$noisy_cdom <- noise + (0.05/1*exp$depth)*exp(0.004*(exp$depth+100))
# exp$noisy_cdom <- noise + (0.05/1*exp$depth)*sin(0.04*exp$depth)

ggplot(exp, aes(x = cdom, y = depth)) + geom_path(colour = "red") + 
  scale_y_reverse() + geom_path(data = obs, aes(x = cdom, y = depth),
                                colour = "black") + 
  geom_path(data = exp, aes(x = noisy_cdom, y = depth), colour = "blue")

calibrange <- exp[which(exp$depth == TopDepth):which(exp$depth == MaxDepth),]

linearMod <- lm(fluo ~ noisy_cdom, data = calibrange)
slope_fdom <- coef(linearMod)[[2]]
C <- coef(linearMod)[[1]]

exp$fluo_cor <- exp$fluo - (slope_fdom*exp$noisy_cdom) - C
exp$fluo_cor <- exp$fluo - (slope_fdom*exp$noisy_cdom) - C


ggplot(exp, aes(x = fluo_cor, y =depth)) + geom_path(colour="red") +
  scale_y_reverse() + 
  geom_path(data = exp, aes(x = fluo, y = depth), colour = "black") + 
  geom_path(data = obs, aes(x = fluo, y = depth), colour = "blue") + 
  geom_vline(xintercept = 0)
  
  

error_out = data.frame(depth = obs$depth, error = abs((obs$fluo-exp$fluo_cor)/obs$fluo)*100)
rmse = sqrt(sum((obs$fluo - exp$fluo_cor)^2)/length(obs$fluo))


