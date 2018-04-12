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

# Extracting data from CHLA-equipped ARGO buoys ------------------

filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

ExtractVar<-function(Var,FloatInfo){
  with(FloatInfo,{
    # This function should return a dataframe for variable with value, qc, iprofile and ilevel
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

#Get each profile from ARGO data -----
profiledf<-ldply(as.list(filename),function(file){
  
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
  chladf <- ExtractVar("CHLA_ADJUSTED",FloatInfo)
  if (all(is.na(chladf)) == TRUE) {
    chladf <- ExtractVar("CHLA",FloatInfo)
  }
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
  data.frame(depth        = chladf$depth,
             juld            = chladf$juld,
             chla           = chladf$value,
             qc              = chladf$qc,
             day             = month.day.year(chladf$juld,c(1,1,1950))$day,
             month           = month.day.year(chladf$juld,c(1,1,1950))$month,
             year            = month.day.year(chladf$juld,c(1,1,1950))$year,
             DOY             = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon             = chladf$lon,
             lat             = chladf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
})

#Save
profiledf_raw <- profiledf

#Remove bad data (see most recent QC definition for RT mode) ====> NOTE que l'on peut faire
# test[test$qc == 4,]
setDT(profiledf, keep.rownames = TRUE)[]
bad_data <- as.numeric(profiledf[profiledf$qc == 4,]$rn)
profiledf$chla[bad_data] <- NA

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
  
  ### Direct use of adjusted values if available ============> D'après Antoine, pas le faire pour la calibration totale
  chladf <- ExtractVar("CHLA_ADJUSTED",FloatInfo)
  if (all(is.na(chladf)) == TRUE) {
    chladf <- ExtractVar("CHLA",FloatInfo)
  }
  
  #Chla subset df et final df
  subchladf <- subset(chladf,select=c("depth","juld","value","qc","lon","lat"))
  colnames(subchladf)[which(colnames(subchladf)=="value")]<-"CHLA"
  
  chladf <- ddply(subchladf,~juld,summarize,
                  qc = qc[which.max(CHLA)],
                  depthmax = depth[which.max(CHLA)],
                  maxvalue = CHLA[which.max(CHLA)],
                  integration = sum(CHLA),
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
crit1 <- argodf[(argodf$bottom < 100),]#Remove profiles too shallow 
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

# #REMOVE 2013 ET 2018 (because not FULL year)
# argodf <- argodf[!argodf$year == 2013,]
# argodf <- argodf[!argodf$year == 2018,]

for (i in 1:length(criteriadf$id)){
  profiledf <- profiledf[!(profiledf$id == criteriadf$id[i]),]
}

# #remove profiles from 2013 and 2018
# profiledf <- profiledf[!profiledf$year == 2013,]
# profiledf <- profiledf[!profiledf$year == 2018,]

##new id before gaussian elimination
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

##### REMOVE MONTHS 10 to 1 ####

# test <- profiledf
# test <- test[!test$month == 10,]
# test <- test[!test$month == 11,]
# test <- test[!test$month == 12,]
# test <- test[!test$month == 1,]
# profiledf <- test
# profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

##### Chloro QC test #################
#////////!!!!!!!!!!!!!!!!!!\\\\\\\\\
#Rajouter un NPQ test? See Xing 2012.

### VISU ###
i <- 166
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
plot(y~x, data=tab, type="l")
i <- i + 1

#Graphe reconstruction
i <- 166
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
#tab <- data.frame(x=tmp$depth,y=tmp$chla)
plot(y~x, data=tab, type="l", main = i)
tmp2 <- approx(tab$x,tab$y, method="linear")
#abline(v=seq(0,900,100) , col="grey" , lwd=0.6)
lines(y~x, data=tmp2, type="l", col="red")
i <- i+1


#HERE PROFILES SHOULD BE OK BECAUSE BAD DATA HAVE BEEN REMOVED#

#### Mignot 2011 ####

fgauss <- function(x, Fsurf, Zdemi, Fmax, Zmax, dz){
  Fsurf*exp((-log(2)/Zdemi)*x) + Fmax*exp(-(x-Zmax)^2/dz^2)
}

fsigmoid <- function(x, Fsurf, Zdemi, s){
  Fsurf*(1/(1+exp((Zdemi-x)*s)))
}

# profilechiant <- tmp$juld[1]#id = 648
# save(profilechiant,file="profilechiant.Rda")
load("profilechiant.Rda")

gaussiandf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  #gaussiandf <- ldply(as.list(1:400), function(i){
  #i<-620
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  #depthindex <- which.min(tmp$depth <= 80)
  depthindex <- which.min(tmp$depth <= 100)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  maxindex <- which.max(tmp$chla)
  #Parameters estimation
  nonNAindex <- which(!is.na(tmp$chla))
  firstnonNA <- nonNAindex[1]
  Fsurf <- tmp$chla[firstnonNA]
  Zmax<- tmp$depth[maxindex]
  Fmax <- tmp$chla[maxindex]
  Zdemi <- tmp$depth[which(tmp$chla <= Fsurf/2)[1]]
  if(is.na(Zdemi)){#cas où chloro surface est déjà faible et que l'on n'atteint pas Fsurf/2 sur les premiers 100m
    #c'est du code 'à la bourrain' mais bon.. nls2 ne veut pas faire son travail
    indexZ <- which(tmp$chla < Fsurf)
    Zdemi <- tmp$depth[indexZ[which.max(indexZ > maxindex)]]
  }
  dz <- 10
  
  if(tmp$juld[1] == profilechiant){
    dz <- 15
  }
  
  res <- nlsLM(y ~ fgauss(x, Fsurf, Zdemi, Fmax, Zmax, dz),
               start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
               data=tab, control = nls.control(maxiter=1000,
                                               minFactor = 1/2048, warnOnly=T))
  
  
  # i <- i+1
  # # v <- summary(res)$parameters[,"Estimate"]
  # # 
  data.frame(Fsurf = coef(res)["Fsurf"], Zdemi = coef(res)["Zdemi"],
             Fmax = coef(res)["Fmax"], Zmax = coef(res)["Zmax"], dz = coef(res)["dz"], id=i,
             juld=tmp$juld[1], file=tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
  # # 
  # # 
  # depthindex <- which.min(tmp$depth <= 100)
  # tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  # plot(y~x, data=tab, type="l")
  # x <- tab$x
  # 
  # 
  # lines(x = tab$x, fgauss(x,coef(res)["Fsurf"], coef(res)["Zdemi"],
  #                         coef(res)["Fmax"], coef(res)["Zmax"], coef(res)["dz"]),
  #                           col=5, add=T, xlim=range(tab$x))
  # 
  # lines(x = tab$x, fsigmoid(x,coef(res)["Fsurf"], coef(res)["Zdemi"],
  #                    coef(res)["s"]),col=3, add=T, xlim=range(tab$x))
  
  
})

sigmoidf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  # sigmoidf <- ldply(as.list(250:600), function(i){
  # i<-164
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  #depthindex <- which.min(tmp$depth <= 80)
  depthindex <- which.min(tmp$depth <= 100)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  maxindex <- which.max(tmp$chla)
  #Parameters estimation
  nonNAindex <- which(!is.na(tmp$chla))
  firstnonNA <- nonNAindex[1]
  Fsurf <- tmp$chla[firstnonNA]
  Zdemi <- tmp$depth[which(tmp$chla <= Fsurf/2)[1]]
  if(is.na(Zdemi)){#cas où chloro surface est déjà faible et que l'on n'atteint pas Fsurf/2 sur les premiers 100m
    #c'est du code 'à la bourrain' mais bon.. nls2 ne veut pas faire son travail
    indexZ <- which(tmp$chla < Fsurf)
    Zdemi <- tmp$depth[indexZ[which.max(indexZ > maxindex)]]
  }
  #s <- (-Fsurf/2)/(Zdemi - tmp$depth[firstnonNA])
  i1 <- which.min(tmp$depth <= Zdemi - 5)
  i2 <- which.min(tmp$depth <= Zdemi + 5)
  d1 <- tmp$depth[i1]
  d2 <- tmp$depth[i2]
  s <- -(Fsurf/2)/(d2-d1)
  
  if (!(Fsurf > 3*tmp$chla[length(tmp$chla)])){
    s <- -2.5
  }
  # st1 <- expand.grid(Fsurf = Fsurf,
  #                    Zdemi= Zdemi,
  #                    s = c(-1,-0.5,-0.1,-0.01,-0.005,-0.001))
  # 
  # st1 <- expand.grid(Fsurf = Fsurf,
  #                    Zdemi= Zdemi,
  #                    s = c(-1,-0.5))
  # 
  # st1 <- c(Fsurf, Zdemi, -1)
  # 
  # mod <- nls2(y ~  fsigmoid(x, Fsurf, Zdemi, s),
  #             start = st1, data=tab, nls.control(maxiter=1000,
  #             minFactor = 1/2048, warnOnly=T), na.action=na.omit,
  #             algorithm = "brute-force", trace=T)
  # 
  # 
  # fo <- y ~ Fsurf*(1/(1+exp((Zdemi-x)*s)))
  # 
  # st1 <- expand.grid(Fsurf = Fsurf,
  #                    Zdemi= Zdemi,
  #                    s = c(-1,-0.5,-0.1,-0.01,-0.005,-2))
  # 
  # fm <- nls2(fo, start = st1,  alg = "plinear-random")
  # fm <- nls2(fo, start = st1,  alg = "random")
  
  # 
  # mod <- nls2(y ~  fsigmoid(x, Fsurf, Zdemi, s),
  #             start = st1,
  #             algorithm = "brute-force")
  # 
  # fo <- y ~ Const + B * (x ^ A)
  # # pass our own set of starting values
  # # returning result of brute force search as nls object
  # st1 <- expand.grid(Const = seq(-100, 100, len = 4),
  #                    B = seq(-100, 100, len = 4), A = seq(-1, 1, len = 4))
  # mod1 <- nls2(fo, start = st1, algorithm = "brute-force")
  
  res <- nlsLM(y ~ fsigmoid(x, Fsurf, Zdemi, s),
               start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
               data=tab, control = nls.control(maxiter=1000,
                                               minFactor = 1/2048, warnOnly=T))
  
  
  #i <- i+1
  # # v <- summary(res)$parameters[,"Estimate"]
  # # 
  data.frame(Fsurf = coef(res)["Fsurf"], Zdemi = coef(res)["Zdemi"], s=coef(res)["s"], id=i,
             juld=tmp$juld[1], file=tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
  
  # depthindex <- which.min(tmp$depth <= 100)
  # tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  # plot(y~x, data=tab, type="l")
  # x <- tab$x
  # lines(x = tab$x, fsigmoid(x,coef(res)["Fsurf"], coef(res)["Zdemi"],
  #                    coef(res)["s"]),col=3, add=T, xlim=range(tab$x))
  
})

#visualisation triplot
i <- 162
tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
                          sigmoidf$s[i]),col=2, add=T, xlim=range(tab$x), lwd = 2)
lines(x = tab$x, fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
                        gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i + 1

# #version en log base 10
# Rcoefdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
#   tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
#   depthindex <- which.min(tmp$depth <= 100)
#   tab <- data.frame(x=tmp$depth[1:depthindex],y=log10(tmp$chla[1:depthindex]))
#   x <- tab$x
#   logsigmoid <- log10(fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
#                                sigmoidf$s[i]))
#   loggauss <- log10(fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
#                            gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i]))
#   # plot(y~x, data=tab, type="l", lwd = 2)
#   # lines(x = tab$x, logsigmoid, col=2, add=T, xlim=range(tab$x), lwd = 2)
#   # lines(x = tab$x, loggauss, col=4, add=T, xlim=range(tab$x), lwd = 2)
#   
#   mean_data <- mean(tab$y, na.rm=T)
#   ss_tot <- sum((tab$y - mean_data)^2, na.rm=T)
#   ss_res_sigmoid <- sum((logsigmoid - tab$y)^2, na.rm = T)
#   ss_res_gaussian <- sum((loggauss - tab$y)^2, na.rm = T)
#   rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)
#   rcoef_gaussian <- 1-(ss_res_gaussian/ss_tot)
#   data.frame(rcoef_sigmoid = rcoef_sigmoid, rcoef_gaussian = rcoef_gaussian)
# })

#version NON log 10
RcoefdfV2 <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  #i <- 1
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  depthindex <- which.min(tmp$depth <= 100)#on peut modifier cette profondeur mais alors il faut regarder pour les fit functions prises au-dessus.... (et avec 100 ça marche)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  x <- tab$x
  sigmoid <- fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
                      sigmoidf$s[i])
  gauss <- fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
                  gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i])
  # plot(y~x, data=tab, type="l", lwd = 2)
  # lines(x = tab$x, sigmoid, col=2, add=T, xlim=range(tab$x), lwd = 2)
  # lines(x = tab$x, gauss, col=4, add=T, xlim=range(tab$x), lwd = 2)
  
  mean_data <- mean(tab$y, na.rm=T)
  ss_tot <- sum((tab$y - mean_data)^2, na.rm=T)
  ss_res_sigmoid <- sum((sigmoid - tab$y)^2, na.rm = T)
  ss_res_gaussian <- sum((gauss - tab$y)^2, na.rm = T)
  rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)
  rcoef_gaussian <- 1-(ss_res_gaussian/ss_tot)
  #i <- i +1
  data.frame(rcoef_sigmoid = rcoef_sigmoid, rcoef_gaussian = rcoef_gaussian, id = tmp$id[1], juld = tmp$juld[1],
             file = tmp$Platform[1], lon = tmp$lon[1], lat = tmp$lat[1],
             day = tmp$day[1], month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1])
})

########## BIG QUESTION ##########

#Pour le moment je fais un fit sur les premiers 100 m car je fais l'hypothèse que l'on est en
#bonne approximation à 0 pour la [chloro] jusqu'à 100.
#Après... si on fais le QC regional test et qu'on règle les profils de chloro de la black sea
#on pourrait le faire sur l'entièreté des profils.. MAIS est-ce que cela vaut le coup? 
#ÇA FERAIT DES REPORTS DANS LES DATES, ÇA M'ARRANGE PAS, IL SERAIT PLUS SAGE DE PROUVER QUE C'EST
#UNE BONNE HYPOTHÈSE LES PREMIERS 100M.

##################################

#classification des profiles

R_crit <- 0.8 #80% de variance expliquée (Navarro c'était 90%)

classifdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  
  #tmp <- profiledf[profiledf$id==i,]
  tmp <- gaussiandf[gaussiandf$id==i,]
  
  if(RcoefdfV2$rcoef_sigmoid[i] & RcoefdfV2$rcoef_gaussian[i] < R_crit){
    classif <- "other"
    score <- "NA"#on s'en fout de leur score
  } 
  else if(RcoefdfV2$rcoef_sigmoid[i] >= RcoefdfV2$rcoef_gaussian[i]){
    classif <- "sigmoid"
    score <- RcoefdfV2$rcoef_sigmoid[i]
  }
  else{
    classif <- "gaussian"
    score <- RcoefdfV2$rcoef_gaussian[i]
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

#Navarro critère 1 is back ! Remove DCM (from gaussian designated profiles) below 1m
#actually, I think it would be interesting to remove DCM peut-être encore < 5m..
#purement subjectif

gaussprofiles <- transform(gaussprofiles,id=as.numeric(factor(juld)))

#établir un autre critère pour virer les profils qui n'ont pas l'allure d'un DCM
#du style : la valeur max de chloro se trouve au début du profil etc...
# TROUVER UN TRUC POUR DISCRIMINER CES PROFILS QUI ONT PASSÉ LE TEST DE LA SIGMOIDE 
# MALENCONTREUSEMENT
testcount <- length(which(gaussprofiles$score < 0.975))#pas forcément le top mais bon ça permet
# d'éliminer encore des profils que ne ressemblent pas trop à une gaussienne mais bien plus à une
#sigmoïde

#get gaussian 'only' profiles
gaussianprofdf <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- profiledf[profiledf$juld == gaussprofiles$juld[i],]
  data.frame(depth = tmp$depth, juld = tmp$juld, chla = tmp$chla, qc = tmp$qc,
             day = tmp$day, month = tmp$month, year = tmp$year, DOY = tmp$DOY,
             Platform = tmp$Platform, lon = tmp$lon, lat = tmp$lat)
})

gaussianprofdf <- transform(gaussianprofdf,id=as.numeric(factor(juld)))

gaussdata <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- gaussiandf[gaussiandf$juld == gaussprofiles$juld[i],]
  data.frame(Fsurf = tmp$Fsurf, Zdemi = tmp$Zdemi, Fmax = tmp$Fmax,
             Zmax = tmp$Zmax, dz = tmp$dz, juld = tmp$juld, file = tmp$file,
             lon = tmp$lon, lat = tmp$lat, day = tmp$day, month = tmp$month,
             year = tmp$year, DOY = tmp$DOY)
})


gaussdata <- transform(gaussdata,id=as.numeric(factor(juld)))


# GAUSSIAN VISU ONLY
i <-100
tmp <- profiledf[profiledf$juld == gaussprofiles$juld[i],]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
x <- tab$x
plot(y~x, data=tab, type="l", lwd = 2, main=gaussprofiles$score[i], sub = 'test')
lines(x = tab$x, fgauss(x,gaussdata$Fsurf[i], gaussdata$Zdemi[i],
                        gaussdata$Fmax[i], gaussdata$Zmax[i], gaussdata$dz[i]),
      col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i +1

#fusion density et depth R codes -> need the find the density associated to the depth DCM....

#reorder and sort by filename
gaussdata <- gaussdata[order(gaussdata$file),]

#function (code pourri -> bcp de problèmes bizarres.......... /-|-\"  
rho_dcm_from_depth_dcm <- function(file, borne_inf, borne_sup){
    
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

  psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    psaldf <- ExtractVar("PSAL",FloatInfo)
  }
  
  df <- ldply(as.list(borne_inf:borne_sup), function(i){

  depth_dcm <- gaussdata$Zmax[i]
  julian_day <- gaussdata$juld[i]

  psaltmp <- psaldf[psaldf$juld == julian_day,]
  temptmp<- tempdf[tempdf$juld == julian_day,]

  #modifier la ligne suivante pour prendre la valeur la plus proche !!!! 
  #indexdepth <- which(temptmp$depth >= depth_dcm)[1]#note que ce n'est pas la meilleure des méthodes vu les différentes vertical schemes des argos
  indexdepth <- which.min(abs(temptmp$depth - depth_dcm))
  psaltmp <- psaltmp[indexdepth,]
  temptmp <- temptmp[indexdepth,]

  psal <- gsw_SA_from_SP(psaltmp$value,depth_dcm,psaltmp$lon,psaltmp$lat)
  temp <- gsw_CT_from_t(psaltmp$value,temptmp$value,depth_dcm)
  rho_anomaly <- gsw_sigma0(psal,temp)
  #i <- i + 1
  data.frame(rho_dcm = rho_anomaly, filename = gaussdata$file[i],
             lon = gaussdata$lon[i], lat = gaussdata$lat[i],
             day = tmp$day[i], month = tmp$month[i], year = tmp$year[i],
             DOY = tmp$DOY[i], juld = tmp$juld[i])
  })
}

index_dcm <- which(duplicated(gaussdata$file) == FALSE)

rho_dcm1 <- rho_dcm_from_depth_dcm(gaussdata$file[index_dcm[1]], index_dcm[1], index_dcm[2]-1)
rho_dcm2 <- rho_dcm_from_depth_dcm(gaussdata$file[index_dcm[2]], index_dcm[2], index_dcm[3]-1)
rho_dcm3 <- rho_dcm_from_depth_dcm(gaussdata$file[index_dcm[3]], index_dcm[3], index_dcm[4]-1)
rho_dcm4 <- rho_dcm_from_depth_dcm(gaussdata$file[index_dcm[4]], index_dcm[4], length(gaussdata$juld))
rho_dcm <- rbind(rho_dcm1, rho_dcm2, rho_dcm3, rho_dcm4)

#visualisation
i <- 414
tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
x <- tab$x
sigmoid <- fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
                    sigmoidf$s[i])
gauss <- fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
                gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i])
plot(y~x, data=tab, type="l", lwd = 2)
lines(x = tab$x, sigmoid, col=2, add=T, xlim=range(tab$x), lwd = 2)
lines(x = tab$x, gauss, col=4, add=T, xlim=range(tab$x), lwd = 2)
i <- i+ 1 


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
    i <- i + 1
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
  
  # tempdf <- tempdf[!tempdf$qc == 4,]
  # tempdf <- tempdf[!tempdf$qc == 3,]
  
  tempdf <- var_average(tempdf)
  #instead of removing bad lines, put NA's ---> will impact graphics
  setDT(tempdf, keep.rownames = TRUE)[]
  bad_data <- as.numeric(tempdf[tempdf$qc == 4,]$rn)
  tempdf$value[bad_data] <- NA
  bad_data <- as.numeric(tempdf[tempdf$qc == 3,]$rn)
  tempdf$value[bad_data] <- NA
  
  
  psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    psaldf <- ExtractVar("PSAL",FloatInfo)
  }
  
  # psaldf <- psaldf[!psaldf$qc == 4,]
  # psaldf <- psaldf[!psaldf$qc == 3,]
  
  psaldf <- var_average(psaldf)
  setDT(psaldf, keep.rownames = TRUE)[]
  bad_data <- as.numeric(psaldf[psaldf$qc == 4,]$rn)
  psaldf$value[bad_data] <- NA
  bad_data <- as.numeric(psaldf[psaldf$qc == 3,]$rn)
  psaldf$value[bad_data] <- NA
  
  #Extract matching rows temperature and salinity data frames
  tempdf <- match_df(tempdf, psaldf, on = c("depth", "juld", "aprofile", "alevel"))
  psaldf <- match_df(psaldf, tempdf, on = c("depth", "juld", "aprofile", "alevel"))
  tempdf <- match_df(tempdf, chladf, on = c("depth", "juld"))
  psaldf <- match_df(psaldf, chladf, on = c("depth", "juld"))
  chladf <- match_df(chladf, psaldf, on = c("depth", "juld"))#135 lignes de plus que ce à quoi je m'attends...
  #24101 pour tempdf et psaldf vs 24236 pour chladf et finaldf (voir après)

  #Sub data frames
  subtempdf <- subset(tempdf, select = c("value","depth","aprofile", "alevel","juld","lon","lat"))
  colnames(subtempdf)[which(colnames(subtempdf)=="value")]<-"TEMP"
  subpsaldf <- subset(psaldf, select = c("value","depth","aprofile","juld"))
  colnames(subpsaldf)[which(colnames(subpsaldf)=="value")]<-"PSAL"
  subchladf <- subset(chladf, select = c("value","depth","aprofile", "alevel","juld","qc"))
  colnames(subchladf)[which(colnames(subchladf)=="value")]<-"CHLA"
  colnames(subchladf)[which(colnames(subchladf)=="aprofile")]<-"profileCHLA"
  colnames(subchladf)[which(colnames(subchladf)=="alevel")]<-"levelCHLA"
  joindf <- join(subtempdf,subpsaldf, by = c("depth", "aprofile", "juld"))
  colnames(joindf)[which(colnames(joindf)=="aprofile")]<-"profileTEMP/PSAL"
  colnames(joindf)[which(colnames(joindf)=="alevel")]<-"levelTEMP/PSAL"
  subchladf <- match_df(subchladf,joindf, on = c("depth", "juld"))
  finaldf <- join(subchladf,joindf, by = c("depth","juld"))

  #use of gsw package
  #NEED CONVERSION OF PRACTICAL SALINITY TO ABSOLUTE SALINITY BEFORE USING THE FOLLOWING FUNCTION
  #NEED CONVERSION OF IN-SITU TEMPERATURE TO CONSERVATIVE TEMPERATURE
  psal <- gsw_SA_from_SP(finaldf$PSAL,finaldf$depth,finaldf$lon,finaldf$lat)
  temp <- gsw_CT_from_t(psal,finaldf$TEMP,finaldf$depth)
  rho_anomaly <- gsw_sigma0(psal,temp)
  
  rhodf <- data.frame(Depth= finaldf$depth,
                      rho_anomaly = rho_anomaly,  
                      CHLA = finaldf$CHLA,
                      TEMP = temp,
                      PSAL = psal,
                      juld            = finaldf$juld,
                      day             = month.day.year(finaldf$juld,c(1,1,1950))$day,
                      month           = month.day.year(finaldf$juld,c(1,1,1950))$month,
                      year            = month.day.year(finaldf$juld,c(1,1,1950))$year,
                      lon             = finaldf$lon,
                      lat             = finaldf$lat)
}

densityprofiledf<-ldply(as.list(filename),function(file){
  
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
  
  FloatInfo<-list(N_PROF=N_PROF,
                  N_LEVELS=N_LEVELS,
                  juld = juld,
                  pres=pres,
                  lon=lon,
                  lat=lat,
                  cycle_number = cycle_number
  )
  
  ### Direct use of adjusted values if available
  chladf <- ExtractVar("CHLA_ADJUSTED",FloatInfo)
  if (all(is.na(chladf)) == TRUE) {
    chladf <- ExtractVar("CHLA",FloatInfo)
  }
  
  #chladf <- chladf[!chladf$qc == 4,]
  
  #instead of removing bad lines, put NA's ---> will impact graphics
  setDT(chladf, keep.rownames = TRUE)[]
  bad_data <- as.numeric(chladf[chladf$qc == 4,]$rn)
  chladf$value[bad_data] <- NA
  
  #Extraction de l'anomalie de densité potentielle
  densitydf <- ExtractDensity(chladf, FloatInfo)
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
  data.frame(depth         =densitydf$Depth,
            density       = densitydf$rho_anomaly,
             juld            = densitydf$juld,
             chla           = densitydf$CHLA,
             day             = month.day.year(densitydf$juld,c(1,1,1950))$day,
             month           = month.day.year(densitydf$juld,c(1,1,1950))$month,
             year            = month.day.year(densitydf$juld,c(1,1,1950))$year,
             DOY             = as.integer(strftime(as.Date(densitydf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon             = densitydf$lon,
             lat             = densitydf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
})

densityprofiledf_raw <- densityprofiledf

#pour supprimer les doublons qui apparaissent
densityprofiledf <- unique(densityprofiledf)
#Assign a profile number according to the (unique due to two decimals) julian day 
densityprofiledf <- transform(densityprofiledf,id=as.numeric(factor(juld)))

# obsolete now?
# #éventuellement -> introduit un peu d'incertitude mais pas tant que cela
# density_duplicate <- which(duplicated(densityprofiledf$density))
# densityprofiledf <- densityprofiledf[-density_duplicate,]

#get gaussian 'only' profiles
densitygaussianprofilesdf <- ldply(as.list(1:length(gaussprofiles$juld)), function(i){
  tmp <- densityprofiledf[densityprofiledf$juld == gaussprofiles$juld[i],]
  data.frame(depth = tmp$depth, density = tmp$density, juld = tmp$juld, chla = tmp$chla, 
             day = tmp$day, month = tmp$month, year = tmp$year, DOY = tmp$DOY,
             Platform = tmp$Platform)
})

#Density graphics -------
densitygaussianprofilesdf2<-ddply(densitygaussianprofilesdf,~month, transform, season=1*(month %in% c(12,1,2))+
                    2*(month %in% c(3,4,5 ))+
                    3*(month %in% c(6,7,8 ))+
                    4*(month %in% c(9,10,11 )))

#Temporal evolution for 2017
tmp <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2017,]
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

#same as DepthDataExtraction
ggplot(perioddf, aes(x=density, y=chla, group=juld)) +
  geom_line() + facet_grid(~group) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) + coord_flip() + scale_x_reverse()

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
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_viridis() +  scale_x_reverse()

ggplot(densitygaussianprofilesdf, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_gradientn(colours = rainbow(5)) + scale_x_reverse()

#for continuous scale
ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_viridis() + scale_x_reverse()

#discrete case
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 1] <- "Winter"
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 2] <- "Spring"
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 3] <- "Summer"
densitygaussianprofilesdf2$season[densitygaussianprofilesdf2$season == 4] <- "Autumn"

ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
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
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_x_reverse()

#2017

tmp <- densitygaussianprofilesdf2[densitygaussianprofilesdf2$year == 2017,]

ggplot(tmp, aes(x=density, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_x_reverse()

ggplot(densitygaussianprofilesdf2, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
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
  geom_line() + facet_grid(~year) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) + ylim(0,4) +
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
  xlim(80, 0) +ylim(0,4) +
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
  xlim(80, 0) +ylim(0,4) +
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
  xlim(80, 0) +ylim(0,4) +
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

#Some statistics on data for DIVAND --------------