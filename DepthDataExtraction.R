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

for (i in 1:length(criteriadf$id)){
  profiledf <- profiledf[!(profiledf$id == criteriadf$id[i]),]
}

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

############
# TEST PAR #
############

pardf <- ExtractVar("DOWNWELLING_PAR_ADJUSTED",FloatInfo)
if (all(is.na(pardf)) == TRUE) {
  pardf <- ExtractVar("DOWNWELLING_PAR",FloatInfo)
}

pardf <- transform(pardf,id=as.numeric(factor(juld)))

#Calcul des isolumes 10 et 1%

isolumedf <- ldply(as.list(1:length(unique(pardf$id))), function(i){
  tmp <- pardf[pardf$id==i,]
  tmp <- tmp[!(tmp$value < 0),]#remove anomalies
  top <- tmp$value[1]
  iso10 <- top/10
  iso10 <- tmp$depth[which.min(tmp$value >= iso10)]
  iso1 <- top/100
  iso1 <- tmp$depth[which.min(tmp$value >= iso1)]
  data.frame(top = top, iso10 = iso10, iso1= iso1)
})

#####################
# Gaussian Fit TEST #
#####################

# i <- 14
# tmp <- profiledf[profiledf$id == i,]
# 
# plot(x=tmp$chla,y=tmp$depth, ylim = rev(c(0,80)))
# plot(x=b$chla,y=b$depth, ylim = rev(c(0,80)), type="l")
# 
# ggplot(tmp, aes(x=depth, y=chla)) +
#   geom_line() +
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(80, 0) + #ylim(0,4) +
#   coord_flip() 
# 
# d <- SMA(tmp$chla,5)
# e <- SMA(tmp$depth,5)
# plot(x=d,y=tmp$depth, ylim = rev(c(0,80)), type="l")

#NLS + gauss from Navarro
gaussdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  #i<-47
  tmp <- profiledf[profiledf$id==i,]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  maxindex <- which.max(tmp$chla)
  depthmax <- tmp$depth[maxindex]
  bg <- tmp$chla[filterpoints]/2#first background test -> chla[n] dépend du filtre choisi..
  bgindex <- which.max(tmp$chla <= bg)
  h <- sum(tmp$chla[1:bgindex], na.rm=T)
  res <- nls( y ~ b0 + h/sigma*sqrt(2*pi)*exp(-1/2*(x-depthmax)^2/sigma^2),
               start = c(h=h, sigma = 5, b0=bg), data=tab, 
              nls.control(maxiter=1000, minFactor = 1/1024,
                          warnOnly=T), na.action=na.omit)
  v <- summary(res)$parameters[,"Estimate"]
  data.frame(depthmax = depthmax, b0 = v[3], h = v[1], sigma = v[2], id=i,
             juld=tmp$juld[1])
  #i <- i+1
})

gaussplot <- function(id, gaussdf, profiledf){
  tmp <- profiledf[profiledf$id==i,]
  depthindex <- which.min(tmp$depth <= 80)
  tab <- data.frame(x=tmp$depth[1:depthindex],y=tmp$chla[1:depthindex])
  plot(y~x, data=tab, type="l")
  plot(function(x) gaussdf$b0[id] + gaussdf$h[id]/gaussdf$sigma[id]*sqrt(2*pi)*exp(-1/2*(x-gaussdf$depthmax[id])^2/gaussdf$sigma[id]^2),
       col=4, add=T, xlim=range(tab$x))
}

i <- 425
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
  data.frame(depth=x, gaussfit=y, obs=tmp$chla[1:depthindex], id=i)
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
  data.frame(FVU=fvu)
}) 

removeid <- which(FVUdf$FVU > 0.1)

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

#############
  
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
ggplot(nprofiledf, aes(x=depth, y=chla, color=month, group=juld)) +
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
ggplot(nprofiledf, aes(x=depth, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
  xlim(80, 0) +ylim(0,4) +
  coord_flip() + scale_color_viridis()

# ggplot(profiledf, aes(x=-depth, y=chla, color=month, group=juld)) +
#   geom_line() + 
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(-50, -10) + coord_flip() + scale_color_gradientn(colours = rainbow(5))

#seasonal analysis (clearer view)
profiledf2<-ddply(nprofiledf,~month, transform, season=1*(month %in% c(12,1,2))+
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

#########################################################################################"

#Visualization of distribution density for different period (showing different approach)
ggplot(argodf, aes(x = depthmax, fill=factor(month))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))

# function "cut" provide a partioning based on variable month in 4 class
ggplot(argodf, aes(x = depthmax, fill=cut(month,4))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0)) + coord_flip()

# --> Il semble qu'il y ait bien un cycle saisonnier marqué
# Maintenant je voudrais définir les saisons à mon gout et voir si ça correspond d'année en année

# Nasty code to get a user-defined season
argodf2<-ddply(argodf,~month, transform, season=1*(month %in% c(12,1,2))+
                                                2*(month %in% c(3,4,5 ))+
                                                3*(month %in% c(6,7,8 ))+
                                                4*(month %in% c(9,10,11 )))

#pour tenir compte des saisons
# argodf4 <- subset(argodf2, select=c("lon","lat","juld","depthmax","season"))
# 
# argodf4 <- argodf4[!(argodf4$depthmax > 100),]
# #argodf4 <- argodf4[!(argodf4$depthmax < 2),]
# 
# write.table(argodf4, file="depthInput_season.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

argodf5 <- subset(argodf2, select=c("lon","lat","juld","depthmax","DOY"))
argodf5 <- argodf5[!(argodf5$depthmax > 100),]
write.table(argodf5, file="depthInput_DOY.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = TRUE)

argodf2$season[which(argodf2$season == 1)]<-"Winter"
argodf2$season[which(argodf2$season == 2)]<-"Spring"
argodf2$season[which(argodf2$season == 3)]<-"Summer"
argodf2$season[which(argodf2$season == 4)]<-"Autumn"

# Le fait de le définir en "factor" va faire qu'il va respecter l'ordre imposé
argodf2$season<-factor(argodf2$season,levels = c("Winter","Spring","Summer","Autumn"))

ggplot(argodf2, aes(x = depthmax, fill=factor(year))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(.~season,scales = "free")

# --> interessant, pour la bouée 6901866, la distribution en été et fort semblable pour 2015 et 2017 et bien différentes pour 2016...
# --> Dû à la variabilité spatiale ? 

#Seasonal distribution (independent of the year)
ggmap(myMap) + 
  geom_point(aes(x=lon, y=lat, color = season), data = argodf2 , alpha = .4) + 
  facet_grid(~season)

ggplot(argodf2, aes(x = depthmax, fill=factor(year))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(.~cut(lon,2),scales = "free")

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

---------------------------------------------------------------------------------
  # Extracting data from WOD ----------------------------------------------------
---------------------------------------------------------------------------------
  #Measured by CTD and OSD... -> En attendant de résoudre le problème de la base
  #de données EMODnet (807 datasets dispos..) -> OK on garde la WOD??
  
  
  
### ATTENTION : Utilisateur doit mettre ses fichiers ship-based dans un 
### sous-répertoire sous le nom /Merged
  
wd <- getwd()

#Creation of a list of files
tmp <- list.files(path = paste0(wd,"/Merged"), pattern="*.nc")
#tmp <- c("wod_011518869O.nc","wod_011518871O.nc","wod_011520385O.nc","wod_011527542O.nc")

#fonction qui va éliminer certains fichiers qui ne répondent pas à des critères
#que nous avons établis sur certaines bases
filedf <- ldply(as.list(tmp), function(file){
  bin <- 1#1 if criteria are ok, 0 otherwise
  ncfile <- nc_open(paste0(wd,"/Merged/",file), write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  chla <- ncvar_get(ncfile,"Chlorophyll")
  nsamples <- length(chla)
  pressure <- ncvar_get(ncfile,"Pressure")
  maxpressure <- pressure[nsamples]
  lat <- ncvar_get(ncfile, "lat")
  lon <- ncvar_get(ncfile,"lon")
  nc_close(ncfile)
  
  #criteria
  if (nsamples < 5){
    bin <- 0
  }
  else if (maxpressure < 60){#pour faire face au problème relevé par Arthur.
    bin <- 0
  }
  else if (length(chla[chla < 0]) != 0){
    bin <- 0
  }
  #Pour virer le plateau continental nord-ouest (et aussi un peu sur l'extrême ouest)
  else if (lat > 44.7 | lon < 28.4 | (lat > 43.7 & lon < 31)){
    bin <- 0
  }

  #list(filename = file, bin = bin)
  #filedf <- data.frame(filename = file, bin = bin)
  data.frame(filename = file, bin = bin)
  # subfiledf <- subset(filedf, bin == 1)
  # subfiledf <- subset(subfiledf, select = "filename")
  # as.list(subfiledf)
  })
  
subfiledf <- subset(filedf, bin == 1)
subfiledf <- subset(subfiledf, select = "filename")
filelist <- as.list(subfiledf$filename) 

shipdf <- ldply(as.list(filelist),function(file){
  
  ncfile <- nc_open(paste0(wd,"/Merged/",file), write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  chla <- ncvar_get(ncfile,"Chlorophyll")
  depth <- ncvar_get(ncfile,"Pressure")
  lat <- ncvar_get(ncfile,"lat")
  long <- ncvar_get(ncfile,"lon")
  time <- ncvar_get(ncfile,"time")
  id <- ncvar_get(ncfile,"wod_unique_cast")
  nc_close(ncfile)
  data.frame(Chlorophyll = chla, Depth = depth, day = month.day.year(time,c(1,1,1770))$day,
             month = month.day.year(time,c(1,1,1770))$month, year = month.day.year(time,c(1,1,1770))$year,
             Latitude = lat, Longitude = long, Platform = id)
  
})

#Final data frame
shipdff <- ddply(shipdf,~Platform,summarize, 
                 depthmax = Depth[which.max(Chlorophyll)],
                 maxvalue = Chlorophyll[which.max(Chlorophyll)],
                 integratedvalue = sum(Chlorophyll), day = day[which.max(Chlorophyll)],
                 month = month[which.max(Chlorophyll)], year = year[which.max(Chlorophyll)],
                 lon = Longitude[which.max(Chlorophyll)], lat = Latitude[which.max(Chlorophyll)],
                 type = "ship-based")#ship-based car pas d'infos précises sur WOD 

#Ajout des saisons                 
shipdf2<-ddply(shipdff,~month, transform, season=1*(month %in% c(12,1,2))+
                 2*(month %in% c(3,4,5 ))+
                 3*(month %in% c(6,7,8 ))+
                 4*(month %in% c(9,10,11 )))

shipdf2$season[which(shipdf2$season == 1)]<-"Winter"
shipdf2$season[which(shipdf2$season == 2)]<-"Spring"
shipdf2$season[which(shipdf2$season == 3)]<-"Summer"
shipdf2$season[which(shipdf2$season == 4)]<-"Autumn"

shipdf2$season<-factor(shipdf2$season,levels = c("Winter","Spring","Summer","Autumn"))


# #Data distribution
# ggplot(datadff, aes(depthmax)) + geom_density()#note qu'il y a bcp de zéro car il semble y avoir pas mal de sampling en surface... (intérêt?)
# 
# #Spatial Coverage
# myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
# ggmap(myMap) + geom_point(aes(x=Longitude, y=Latitude), data = datadff, alpha = .2)


#Merging of ARGO and ship-based samplings  -------------------------------------------
merged <- merge(argodf2,shipdf2,all=T)
merged <- subset(merged, select = c("depthmax", "maxvalue", "integratedvalue", 
                                    "day", "month", "year", "lon", "lat",
                                    "Platform", "type", "season"))

ggplot(shipdf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(~season,scales = "free")

ggplot(argodf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(~season,scales = "free")

ggplot(shipdf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))

ggplot(argodf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))

# Spatial analysis of filtered ship files -------
#Spatial Coverage
myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
h <- ggmap(myMap) +
geom_point(aes(x=lon, y=lat, color = season), data = merged, alpha = .4) + 
facet_grid(~type)  

#Seasonal Coverage
sc <- ggmap(myMap) + 
  geom_point(aes(x=lon, y=lat, color = season), data = merged, alpha = .4) + 
  facet_grid(type~season) 


library(gganimate)
library(gapminder)
test <- ggmap(myMap) + 
  geom_point(aes(x=lon, y=lat, color = season, frame = year), data = merged, alpha = .5) +
  facet_grid(season~type)

gganimate(test, interval = 3)
#gganimate(test, "general.gif")

#JUST FOR ARGO DATA FOR NOW
#Data export test on argodf2
# argodf3 <- subset(argodf2, select=c("lon","lat","juld","depthmax"))
# 
# #elimination of rows containing NA's (not yet, do bash) or
# # outliers (here defined for depthmax > 100 and depthmax < 5m)
# argodf3 <- argodf3[!(argodf3$depthmax > 100),]
# argodf3 <- argodf3[!(argodf3$depthmax < 2),]#très arbitraire.... -> checker la doc des ARGOS pour le min
# 
# write.table(argodf3, file="depthInput.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

# write.table(argodf3, file="rhoInput.csv", sep=",", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

# write.table(argodf3, file="rhoInput.txt", sep="\t", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

#system("bash depthinputfile.sh")#execution of a shell script to have a 'clean' input file
