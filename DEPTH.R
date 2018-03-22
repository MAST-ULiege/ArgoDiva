library(ncdf4)
#More complex libraries for plotting data
library(ggplot2)
library(plyr)
#for melt function
library(reshape2)
#for julian date manipulation
library(chron)
library(viridis)#same color palette than Julia's default one
library(TTR)#for SMA filter
library(minpack.lm)
library(nls2)
library(nlstools)
#library(dplyr)
#library(signal) -> bof aussi mais  existe
#library(anomalyDetection)#to spot outliers etc --> bof

# Extracting data from ARGO buoys ------------------

#Seuls fichiers qui possèdent des données de CHLORO. Seuls deux d'entrent eux possèdent aussi du CDOM
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")
# filename <- "6901866_Mprof.nc"
# filename <- "7900591_Mprof.nc"
# filename <- "6900807_Mprof.nc"
# filename <- "7900592_Mprof.nc"

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

#Construction of an ARGO-only dataframe
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

  ### Direct use of adjusted values if available
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
                          lat=mean(lat))#mean car la bouée a bougé sur la journée
  
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
             DOY             = as.integer(strftime(as.Date(chladf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             lon             = chladf$lon,
             lat             = chladf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
})


#Assign a profile number according to the (unique due to two decimals) julian day 
argodf <- transform(argodf,id=as.numeric(factor(juld)))
#argodf <- transform(argodf,id=sort(as.numeric(factor(juld))))
  
#############################################################   ALL PROFILES (TFE)
#Construction of an ARGO-only dataframe
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
             #lon             = chladf$lon,
             #lat             = chladf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
})


#Set profiles id
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

#Get profiles id that do not answer to some criteria
crit1 <- argodf[(argodf$depthmax < 1),]#Remove profiles for which DCM depth is below 1m
crit2 <- argodf[(argodf$bottom < 100),]#Remove profiles too shallow (à voir si on garde ou pas ->
#fonction de la NLS...)
crit3 <- argodf[(argodf$depthmax == argodf$bottom),]#Bad profiles (all values at the same depth)
crit4 <- argodf[(argodf$depthmax >= 80),] #car nls se fait sur les premiers 80m
#rmq : pas forcément des fonctions mais plusieurs valeurs de chla pour une même prof... 
# --> unicité? (unique function)
#crit sur profiledf.. -> en lien avec le bond soudain dans les données... -> ok
#il faut faire cela mais comment faire et quelle ligne garder??? -> on peut supprimer les deux
#pour être certain puis on applique un filtre pour avoir une nls qui ne foire plus trop..
criteriadf <- rbind(crit1, crit2, crit3, crit4)
criteriadf <- unique(criteriadf)

#nouveau selection criteria... -> car l'exclusion du dcm < 1m ne peut se faire qu'après le fit..
#la méthode de Navarro est pas super top à ce niveau, la gars n'avait aucune sigmoide ou quoi??

crit2 <- argodf[(argodf$bottom < 100),]#Remove profiles too shallow (à voir si on garde ou pas ->
crit4 <- argodf[(argodf$depthmax >= 100),] #car nls se fait sur les premiers 80m
crit3 <- argodf[(argodf$depthmax == argodf$bottom),]#Bad profiles (all values at the same depth)
criteriadf <- rbind(crit2, crit3, crit4)
criteriadf <- unique(criteriadf)

#### TEST

forRho1 <- as.data.frame(criteriadf$juld)
save(forRho1,file="firstJULDs.Rda")


# firstIDtoremove <- criteriadf$id
# write.table(firstIDtoremove, file="firstIDs.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)
# 
# firstIDtoremove <- criteriadf$juld
# write.table(firstIDtoremove, file="firstJULDs.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

### END TEST

#save
argodf_raw <- argodf

#remove odd profiles (but we keep > 10% FVU gaussian profiles)
for (i in 1:length(criteriadf$id)){
  argodf <- argodf[!(argodf$id == criteriadf$id[i]),]
}

#REMOVE 2013 ET 2018 (because not FULL year)
argodf <- argodf[!argodf$year == 2013,]
argodf <- argodf[!argodf$year == 2018,]


#argodf <- subset(argodf, select=c("lon","lat","juld","depthmax"))
argodf <- subset(argodf, select=c("lon","lat","depthmax"))
#Remove NA's
argodf <- argodf[complete.cases(argodf),]
write.table(argodf, file="ForDivaWEB.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
             row.names = FALSE, col.names = FALSE)


# #TEST HADRIEN
# hadriendf <- argodf
# for (i in 1:length(criteriadf$id)){
#   hadriendf <- hadriendf[!(hadriendf$id == criteriadf$id[i]),]
# }
# 
# hadriendf2 <- subset(hadriendf, select=c("lon","lat","depthmax"))
# write.table(hadriendf2, file="depthInputHADRIEN.csv", sep=",", na = "NA", dec = ".", eol = "\r\n",
#              row.names = FALSE, col.names = TRUE)

###########


for (i in 1:length(criteriadf$id)){
  profiledf <- profiledf[!(profiledf$id == criteriadf$id[i]),]
}

#remove profiles from 2013 and 2018
profiledf <- profiledf[!profiledf$year == 2013,]
profiledf <- profiledf[!profiledf$year == 2018,]


##new id before gaussian elimination
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

#keep non filtered values
rawprofiledf <- profiledf

filterpoints <- 10
# remove depth duplicate
profiledf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]
  depth_duplicate <- which(duplicated(tmp$depth))
  # if(nrow(tmp) != 0){
  #   tmp <- tmp[-depth_duplicate,]
  # }
  if(length(depth_duplicate) != 0){
    tmp <- tmp[-depth_duplicate,]
  }
  meanchla <- SMA(tmp$chla,filterpoints)#Simple Moving Average (remove sharp peaks and smooth the global curve)
  #meandepth <- SMA(tmp$depth,5)
  tmp$chla <- meanchla
  data.frame(tmp = tmp)
})

colnames(profiledf) <- c("depth", "juld","chla","day","month","year","DOY","Platform","type","id")

##### REMOVE MONTHS 10 to 1 ####

test <- profiledf
test <- test[!test$month == 10,]
test <- test[!test$month == 11,]
test <- test[!test$month == 12,]
test <- test[!test$month == 1,]
profiledf <- test
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))



##### Chloro QC test ##### D'ortenzio 2010 (pas seulement)

# Spike and gradient tests -> modifier le threshold? -> Argo QC Chloro 2014 new version 

QCdf <- ldply(as.list(1:length(unique(profiledf$id))), function(j){
  tmp <- profiledf[profiledf$id==j,]
  chla <- tmp$chla
  
  #Range test 
  range <- sapply(1:(length(chla)),function(i){
    if (chla[i] > 50 | chla[i] < -0.1){
      index3 <- i
    }
  })
  
  index3 <- unlist(range)
  
  
  #Spike test (QC test 9) 
  #initialize empty array (test <- ()) si on modifie le code un jour
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
  
  #Gradient test (QC test 11)
  gradient<- sapply(2:(length(chla)-1),function(i){
    Test_Value <- abs(chla[i]-((chla[i+1]+chla[i-1])/2))
    if (Test_Value > 3){  #Treshold_Value <- 3 #3mg/m³ (see d'Ortenzio et al. 2010)
      index2 <- i
    }
  }
)
  index2 <- unlist(gradient)
  
  index <- c(index1, index2, index3, index4)
  index <- unique(index)
  id <- rep(j,length(index))
  data.frame(index = index, id = id)
})

#////////!!!!!!!!!!!!!!!!!!\\\\\\\\\
#Rajouter un NPQ test? See Xing 2012.

i <- 340
i <- 166
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
plot(y~x, data=tab, type="l")
i <- i + 1

#Remove bad data from profiledf
for (i in 1:length(QCdf$id)){
  tmp <- profiledf[profiledf$id == QCdf$id[i],]
  tmp$chla[QCdf$index[i]] <- NA
  profiledf[profiledf$id == QCdf$id[i],] <- tmp
}

#Graphe reconstruction
i <- 56
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 100)
#tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
tab <- data.frame(x=tmp$depth,y=tmp$chla)
plot(y~x, data=tab, type="l", main = i)
tmp2 <- approx(tab$x,tab$y, method="linear")
#abline(v=seq(0,900,100) , col="grey" , lwd=0.6)
lines(y~x, data=tmp2, type="l", col="red")
i <- i+1

tested <- profiledf

# save(tested,file="profiledfaftertests.Rda")
load("profiledfaftertests.Rda")

#repasser le QC une seconde fois (regarder sur i = 166 pour voir -> car 
# des pics se suivent donc ils vont encore rester par après)
#note que c'est un excellent exemple pour montrer l'effet des 2 QC 
# figure à mettre dans le mémoire -> marchera pas pour les pics côte
# à côte car NA pour le calcul des treshold etc donc NOPE une seule fois
# d'ailleurs si plusieurs pics se suivent, why not.

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
  # gaussiandf <- ldply(as.list(1:400), function(i){
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
             juld=tmp$juld[1])
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
             juld=tmp$juld[1])

  # depthindex <- which.min(tmp$depth <= 100)
  # tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  # plot(y~x, data=tab, type="l")
  # x <- tab$x
  # lines(x = tab$x, fsigmoid(x,coef(res)["Fsurf"], coef(res)["Zdemi"],
  #                    coef(res)["s"]),col=3, add=T, xlim=range(tab$x))

})

#visualisation triplot
i <- 100
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

#version en log base 10
Rcoefdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]#gappy data, no interpolation
  depthindex <- which.min(tmp$depth <= 100)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=log10(tmp$chla[1:depthindex]))
  x <- tab$x
  logsigmoid <- log10(fsigmoid(x,sigmoidf$Fsurf[i], sigmoidf$Zdemi[i],
                               sigmoidf$s[i]))
  loggauss <- log10(fgauss(x,gaussiandf$Fsurf[i], gaussiandf$Zdemi[i],
         gaussiandf$Fmax[i], gaussiandf$Zmax[i], gaussiandf$dz[i]))
  # plot(y~x, data=tab, type="l", lwd = 2)
  # lines(x = tab$x, logsigmoid, col=2, add=T, xlim=range(tab$x), lwd = 2)
  # lines(x = tab$x, loggauss, col=4, add=T, xlim=range(tab$x), lwd = 2)
  
  mean_data <- mean(tab$y, na.rm=T)
  ss_tot <- sum((tab$y - mean_data)^2, na.rm=T)
  ss_res_sigmoid <- sum((logsigmoid - tab$y)^2, na.rm = T)
  ss_res_gaussian <- sum((loggauss - tab$y)^2, na.rm = T)
  rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)
  rcoef_gaussian <- 1-(ss_res_gaussian/ss_tot)
  data.frame(rcoef_sigmoid = rcoef_sigmoid, rcoef_gaussian = rcoef_gaussian)
})

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
  data.frame(rcoef_sigmoid = rcoef_sigmoid, rcoef_gaussian = rcoef_gaussian, id = tmp$id[1], juld = tmp$juld[1])
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
  
  data.frame(classif = classif, id = i, juld = tmp$juld[1], depthmax = tmp$Zmax[1], score = score)
})
  
#count profiles
ngauss <- length(which(classifdf$classif == "gaussian"))/length(unique(profiledf$id))*100
nsigmoid <- length(which(classifdf$classif == "sigmoid"))/length(unique(profiledf$id))*100
nother <- length(which(classifdf$classif == "other"))/length(unique(profiledf$id))*100

#gaussian profile (classified) before a new criteria elimination
gaussprofiles <- classifdf[classifdf$classif == "gaussian",]
  
#remove profiles where depthmax (DCM) < 1m
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

 
#visu
i <- 200
tmp <- profiledf[profiledf$juld == gaussprofiles$juld[i],]
depthindex <- which.min(tmp$depth <= 100)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
plot(y~x, data=tab, type="l", lwd = 2, main=gaussprofiles$score[i])
i <- i +1
#fusion density et depth R codes -> need the find the density associated to the depth DCM....


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


#####################
# Gaussian Fit TEST #
#####################

f <- function(x, bg_up, sigma_r, sigma_l, h, zm){
  ifelse(x > zm, h/sigma_r*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_r^2),
         bg_up + (h-bg_up)/sigma_l*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_l^2))
}

gaussTESTdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  i<-50
  tmp <- profiledf[profiledf$id==i,]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  maxindex <- which.max(tmp$chla)
  #Parameters
  zm<- tmp$depth[maxindex]
  #si on garde le filter
  bg_up <- tmp$chla[filterpoints]
  #bg_down = 0
  sigma_r <- 5#r for right
  sigma_l <- 5
  h <- tmp$chla[maxindex]
  res2 <- nls(y ~ f(x, bg_up, sigma_r, sigma_l, h, zm),
             start = c(h=h, sigma_r = 5, sigma_l = 5, bg_up = bg_up,
                       zm = zm), data=tab, nls.control(maxiter=1000,
                  minFactor = 1/2048, warnOnly=T), na.action=na.omit)
  res <- nlsLM(y ~ f(x, bg_up, sigma_r, sigma_l, h, zm),
             start = c(h=h, sigma_r = 5, sigma_l = 5, bg_up = bg_up,
                       zm = zm), data=tab)
  
  v <- summary(res)$parameters[,"Estimate"]
  
  data.frame(depthmax = zm, id=i,
             juld=tmp$juld[1])
  i <- i+1
})

test <- 56
tmp <- profiledf[profiledf$id==test,]
depthindex <- which.min(tmp$depth <= 80)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
maxindex <- which.max(tmp$chla)
#Parameters
zm<- tmp$depth[maxindex]
#si on garde le filter
bg_up <- tmp$chla[filterpoints]
#bg_down = 0
sigma_r <- 5#r for right
sigma_l <- 5
h <- tmp$chla[maxindex]

# f <- function(x, bg_up, sigma_r, sigma_l, h, zm){
#   ifelse(x > zm, h/sigma_r*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_r^2),
#          bg_up + (h-bg_up)/sigma_l*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_l^2))
# }

f <- function(x, bg_up, sigma_r, sigma_l, h, zm){
  ifelse(x > zm, h/sigma_r*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_r^2),
         bg_up + (h-bg_up)/sigma_l*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_l^2))
}

formula <- as.formula(y ~ (x <= zm) * 
      bg_up + (h-bg_up)/sigma_l*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_l^2) + 
  (x > zm) * h/sigma_r*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_r^2))

preview(formula, data = tab$x, start = list(h, bg_up, sigma_l, sigma_r, zm))

res <- nls(y ~ f(x, bg_up, sigma_r, sigma_l, h, zm),
             start = c(h=h, sigma_r = 5, sigma_l = 5, bg_up = bg_up,
            zm = zm), data=tab, nls.control(maxiter=1000,
            minFactor = 1/2048, warnOnly=T), na.action=na.omit)

v <- summary(res)$parameters[,"Estimate"]
depthindex <- which.min(tmp$depth <= 80)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
plot(y~x, data=tab, type="l")
x <- tab$x
lines(x = tab$x, f(x,coef(res)["bg_up"], coef(res)["sigma_r"], coef(res)["sigma_l"],
     coef(res)["h"], coef(res)["zm"]),col=3, add=T, xlim=range(tab$x))


################# TEST VISU #####

f <- function(x, bg_up, sigma_r, sigma_l, h, zm){
  ifelse(x > zm, h/sigma_r*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_r^2),
         bg_up + (h-bg_up)/sigma_l*sqrt(2*pi)*exp(-1/2*(x-zm)^2/sigma_l^2))
}

gaussTESTdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
i<-156
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 80)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
maxindex <- which.max(tmp$chla)
#Parameters
zm<- tmp$depth[maxindex]
#si on garde le filter
bg_up <- tmp$chla[filterpoints]
#bg_down = 0
sigma_r <- 5#r for right
sigma_l <- 5
h <- tmp$chla[maxindex]

nls2(y ~  f(x, bg_up, sigma_r, sigma_l, h, zm),
     start = c(h=h, sigma_r = 10, sigma_l = 10, bg_up = bg_up,
               zm = zm), data=tab, nls.control(maxiter=1000,
                minFactor = 1/2048, warnOnly=T), na.action=na.omit,
     algorithm = "brute-force")

nls2(y ~  f(x, bg_up, sigma_r, sigma_l, h, zm),
     start = , data=tab, nls.control(maxiter=1000,
     minFactor = 1/2048, warnOnly=T), na.action=na.omit,
     algorithm = "brute-force")

mod <- nls2(y ~  f(x, bg_up, sigma_r, sigma_l, h, zm),
     start = st1, data=tab, nls.control(maxiter=500,
     minFactor = 1/2048, warnOnly=T), na.action=na.omit,
     algorithm = "brute-force", trace=T)

st1 <- expand.grid(h = h,
                   sigma_l = seq(1,100, len=20), 
                   sigma_r = seq(1,100, len=20),
                   bg_up = bg_up,
                   zm = zm)

### remove outlier SPIKE ###




######

st2 <- data.frame(h = c(0, 10), sigma_r = c(1, 100), 
                  sigma_l = c(1,100), bg_up = c(0,10), zm = c(1,80))

res2 <- nls(y ~ f(x, bg_up, sigma_r, sigma_l, h, zm),
            start = c(h=h, sigma_r = 9 , sigma_l = 4, bg_up = bg_up,
                      zm = zm), data=tab, nls.control(maxiter=1000,
                                                      minFactor = 1/2048, warnOnly=T), na.action=na.omit)


res <- nlsLM(y ~ f(x, bg_up, sigma_r, sigma_l, h, zm),
             start = c(h=h, sigma_r = 9, sigma_l = 4, bg_up = bg_up,
                       zm = zm), data=tab)

res <- nlsLM(y ~ f(x, bg_up, sigma_r, sigma_l, h, zm),
             start = coef(mod), data=tab)
svd(mod$m$Rmat())[-2]

plot(y~x, data=tab, type="l")
x <- tab$x
lines(x = tab$x, f(x,coef(res)["bg_up"], coef(res)["sigma_r"], coef(res)["sigma_l"],
                   coef(res)["h"], coef(res)["zm"]),col=3, add=T, xlim=range(tab$x))

lines(x = tab$x, f(x,coef(res2)["bg_up"], coef(res2)["sigma_r"], coef(res2)["sigma_l"],
                   coef(res2)["h"], coef(res2)["zm"]),col=3, add=T, xlim=range(tab$x))

i <- i+1
})

########### END test    #####

########## TEST LOGISTIC GROWTH (REVERSE) ######

logisticdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
i<-1
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 80)
tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
tab <- data.frame(x=tmp$depth[filterpoints:depthindex],y=tmp$chla[filterpoints:depthindex])
plot(y~x, data=tab, type="l")
# res2 <- nls(y ~ A/(1+exp(k*(x-d))),
#             start = c(A=1, k = 2, d = 45), data=tab, nls.control(maxiter=1000,
#                                                       minFactor = 1/2048, warnOnly=T), na.action=na.omit)

#nls2(y ~ A/(1+exp(k*(x-d))), data=tab, start = c(A=1, k = 2, d = 45))

res <- nlsLM(y ~ A/(1+exp(k*(x-d))),
             start = c(A=1, k = 2, d = 45), data=tab,
             control = nls.lm.control(maxiter = 500))
#i <- i + 1
preview (y ~ A/(1+exp(k*(x-d))),data=tab, start = c(A=1, k = 2, d = 45),
         variable = 1)
plotfit (res, pch.obs = 1, pch.fit = "+", lty = 1, lwd = 1,
         col.obs = "black", col.fit ="red")
overview(res)
i <- i +1
plot(y~x, data=tab, type="l")
x <- tab$x
lines(x = tab$x, coef(res2)["A"]/(1+exp(coef(res2)["k"]*(x-coef(res2)["d"]))),
      col=4, add=T, xlim=range(tab$x))
lines(x = tab$x, coef(res)["A"]/(1+exp(coef(res)["k"]*(x-coef(res)["d"]))),
      col=2, add=T, xlim=range(tab$x))

v <- summary(res)$parameters[,"Estimate"]

data.frame(A=v[1], k=v[2], d=v[3], id=i,
           juld=tmp$juld[1])
#i <- i+1
})

logvardf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[filterpoints:depthindex],
                    y=tmp$chla[filterpoints:depthindex])
  x <- tab$x
  y <- logisticdf$A[id]/(1+exp(logisticdf$k[id]*(x-logisticdf$d[id])))
  #obs = SMA observations ! and maybe we could also calcule diff avec la
  #gaussienne à partir de la n ème valeur (#filterpoints)
  data.frame(depth=x, fit=y, obs=tmp$chla[filterpoints:depthindex],
             id=i, juld = tmp$juld[1])
  # data.frame(x=x, gaussfit=y[filterpoints:depthindex],
  #            obs=tmp$chla[filterpoints:depthindex], id=i)
})

#marche pas du tout du tout du tout WHHYYYYYY
RMSdf <- ldply(as.list(1:length(unique(logisticdf$id))), function(i){
  tmp <- logvardf[logvardf$id==i,]
  fit <- tmp$fit[filterpoints:length(tmp$fit)]
  obs <- tmp$obs[filterpoints:length(tmp$fit)]
  RMSE <- sqrt(sum((obs-fit)^2)/length(obs))
  data.frame(RMSE=RMSE, juld=tmp$juld[1])
}) 

cordf <- ldply(as.list(1:length(unique(logisticdf$id))), function(i){
  tmp <- logvardf[logvardf$id==i,]
  fit <- tmp$fit[filterpoints:length(tmp$fit)]
  obs <- tmp$obs[filterpoints:length(tmp$fit)]
  cor <- cor(x=fit, y=obs, use="na")
  data.frame(R=cor, juld=tmp$juld[1])
}) 

logisticplot <- function(id, logisticdf, profiledf, RMSdf){
  tmp <- profiledf[profiledf$id==i,]
  time <- tmp$month[filterpoints]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  plot(y~x, data=tab, type="l", main = paste0(time," - ",RMSdf$RMSE[id]))
  #plot(y~x, data=tab, type="l", main = time)
  x <- tab$x
  plot(function(x) logisticdf$A[id]/(1+exp(logisticdf$k[id]*(x-logisticdf$d[id]))),
       col=4, add=T, xlim=range(tab$x))
}

i <- 1
logisticplot(i, logisticdf, profiledf, RMSdf)
i<- i+1

######## END TEST LOGISTIC (REVERSE) ########



#### Test organelli #####

i<-56
tmp <- profiledf[profiledf$id==i,]
depthindex <- which.min(tmp$depth <= 80)
tab <- data.frame(x=tmp$depth[filterpoints:depthindex],y=tmp$chla[filterpoints:depthindex])





#NLS + gauss from Navarro
gaussdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  i<-94
  tmp <- profiledf[profiledf$id==i,]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  maxindex <- which.max(tmp$chla)
  depthmax <- tmp$depth[maxindex]
  zm <- depthmax
  bg <- tmp$chla[filterpoints]/2#first background test -> chla[n] dépend du filtre choisi..
  bgindex <- which.max(tmp$chla <= bg)
  h <- sum(tmp$chla[1:bgindex], na.rm=T)
  res <- nls( y ~ b0 + h/sigma*sqrt(2*pi)*exp(-1/2*(x-depthmax)^2/sigma^2),
               start = c(h=h, sigma = 5, b0=bg), data=tab, 
              nls.control(maxiter=1000, minFactor = 1/2048,
                          warnOnly=T), na.action=na.omit)
  v <- summary(res)$parameters[,"Estimate"]
  data.frame(depthmax = depthmax, b0 = v[3], h = v[1], sigma = v[2], id=i,
             juld=tmp$juld[1])
  #i <- i+1
})

gaussplot <- function(id, gaussdf, profiledf){
  tmp <- profiledf[profiledf$id==i,]
  time <- tmp$month[filterpoints]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  plot(y~x, data=tab, type="l", main = time)
  plot(function(x) gaussdf$b0[id] + gaussdf$h[id]/gaussdf$sigma[id]*sqrt(2*pi)*exp(-1/2*(x-gaussdf$depthmax[id])^2/gaussdf$sigma[id]^2),
     col=4, add=T, xlim=range(tab$x))
  # plot(function(x) v[3] + v[1]/v[2]*sqrt(2*pi)*exp(-1/2*(x-v[4])^2/v[2]^2),
  #      col=4, add=T, xlim=range(tab$x))
}



###############################

i <- 1
gaussplot(i, gaussdf, profiledf)
i<- i+1

gaussvardf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  x <- tab$x
  y <- gaussdf$b0[i] + gaussdf$h[i]/gaussdf$sigma[i]*sqrt(2*pi)*exp(-1/2*(x-gaussdf$depthmax[i])^2/gaussdf$sigma[i]^2)
  #obs = SMA observations ! and maybe we could also calcule diff avec la
  #gaussienne à partir de la n ème valeur (#filterpoints)
  data.frame(depth=x, gaussfit=y, obs=tmp$chla[1:depthindex], id=i, juld = tmp$juld[1])
  # data.frame(x=x, gaussfit=y[filterpoints:depthindex],
  #            obs=tmp$chla[filterpoints:depthindex], id=i)
  })

cordf <- ldply(as.list(1:length(unique(gaussvardf$id))), function(i){
  tmp <- gaussvardf[gaussvardf$id==i,]
  cor <- cor(x=tmp$gaussfit, y=tmp$obs, use="na")
  cov <- cov(x=tmp$gaussfit, y=tmp$obs, use="na")
  data.frame(cor=cor, cov=cov, id=i)
})

#note : je pourrais aussi calculer une corrélation ou autre juste entre
#zm+sigma et zm-sigma 
#keepid <- which(!(cordf$cor < 0.95))

# FRACTION OF VARIANCE UNEXPLAINED
# https://en.wikipedia.org/wiki/Fraction_of_variance_unexplained

FVUdf <- ldply(as.list(1:length(unique(gaussvardf$id))), function(i){
  tmp <- gaussvardf[gaussvardf$id==i,]
  VARobs <- var(x=tmp$obs, na.rm=TRUE)#variance des données que l'on essaie de prédire avec la gaussienne
  gaussfit <- tmp$gaussfit[filterpoints:length(tmp$gaussfit)]
  obs <- tmp$obs[filterpoints:length(tmp$gaussfit)]
  MSE <- sum((obs-gaussfit)^2)/length(obs)
  fvu <- MSE/VARobs
  data.frame(FVU=fvu, juld=tmp$juld[1])
}) 

FVUdf <- transform(FVUdf,id=as.numeric(factor(juld)))
removeid <- which(FVUdf$FVU > 0.1)
forRho2 <- FVUdf$juld[removeid]
forRho1 <- as.data.frame(forRho1)
colnames(forRho1) <- "juld"
forRho2 <- as.data.frame(forRho2)
colnames(forRho2) <- "juld"
test <- rbind(forRho1,forRho2)
save(test,file="data.Rda")
# secondIDtoremove <- FVUdf$id[removeid]
# write.table(secondIDtoremove, file="secondIDs.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

removejuld <- FVUdf$juld[removeid]
for (i in 1:length(removejuld)){
  argodf <- argodf[!argodf$juld == removejuld[i],]
}

argodf <- subset(argodf, select=c("lon","lat","juld","depthmax"))
argodf <- subset(argodf, select=c("lon","lat","depthmax"))
argodf <- argodf[complete.cases(argodf),]
write.table(argodf, file="ForDivaWEB_421.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = FALSE)

#"clean" profiledf -> n for new
nprofiledf <- rawprofiledf#if we want raw values
nprofiledf <- profiledf#if we want smoothed values -> easier interpretation..
for (i in 1:length(removeid)){
  nprofiledf <- nprofiledf[!(nprofiledf$id == removeid[i]),]
}

nprofiledf <- transform(nprofiledf,id=as.numeric(factor(juld)))

#artifice pour float precision problem for density profiles exclusion
# dat <- as.data.frame(unique(nprofiledf$juld))*1000000
# dat <- round(dat, digits= 0)
# write.table(dat, file="dat_juld.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = FALSE)

dat <- subset(nprofiledf, select=c("DOY","year","Platform"))
dat <- unique(dat)
write.table(dat, file="dat_juld.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = TRUE)#note j'ai perdu 1 ligne de données avec cela

#"clean" gaussdf pour voir les courbes
ngaussdf <- gaussdf
for (i in 1:length(removeid)){
  ngaussdf <- ngaussdf[!(ngaussdf$id == removeid[i]),]
}

#resize id
ngaussdf <- transform(ngaussdf,id=as.numeric(factor(juld)))



i <- 434
gaussplot(i, ngaussdf, nprofiledf)
i<- i+1

