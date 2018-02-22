library(ncdf4)
#More complex libraries for plotting data
library(ggplot2)
library(plyr)
#for melt function
library(reshape2)
#for julian date manipulation
library(chron)
#for computing potential density anomalies
library(gsw)
library(viridis)#same color palette than Julia's default one
library(TTR)#for SMA filter

# Extracting data from ARGO buoys ------------------
#potential density anomaly as y axis

#Seuls fichiers qui possèdent des données de CHLORO. Seuls deux d'entrent eux possèdent aussi du CDOM
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")
# filename <- "7900591_Mprof.nc"
# filename <- "6901866_Mprof.nc"

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

psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
if (all(is.na(psaldf)) == TRUE) {
  psaldf <- ExtractVar("PSAL",FloatInfo)
}

#Extract matching rows of a data frame
tempdf <- match_df(tempdf,chladf, on = c("depth", "juld"))
psaldf <- match_df(psaldf,chladf, on = c("depth", "juld"))

#Sub data frames
subtempdf <- subset(tempdf, select = c("value","depth","aprofile", "alevel","juld","lon","lat"))
colnames(subtempdf)[which(colnames(subtempdf)=="value")]<-"TEMP"
subpsaldf <- subset(psaldf, select = c("value","depth","aprofile","juld"))
colnames(subpsaldf)[which(colnames(subpsaldf)=="value")]<-"PSAL"
subchladf <- subset(chladf, select = c("value","depth","aprofile", "alevel","juld"))
colnames(subchladf)[which(colnames(subchladf)=="value")]<-"CHLA"
colnames(subchladf)[which(colnames(subchladf)=="aprofile")]<-"profileCHLA"
colnames(subchladf)[which(colnames(subchladf)=="alevel")]<-"levelCHLA"
joindf <- join(subtempdf,subpsaldf, by = c("depth", "aprofile", "juld"))
colnames(joindf)[which(colnames(joindf)=="aprofile")]<-"profileTEMP/PSAL"
colnames(joindf)[which(colnames(joindf)=="alevel")]<-"levelTEMP/PSAL"
subchladf <- match_df(subchladf,joindf, on = c("depth", "juld"))
finaldf <- join(subchladf,joindf, by = c("depth","juld"))

# #Ascending or descending profile? Or both?
# direction <- ncvar_get(ncfile,"DIRECTION")
# direction <- unlist(strsplit(direction,split=""))
# direction[direction == "A"] <- 1
# direction[direction == "D"] <- 2
# direction <- as.numeric(direction)


#use of gsw package
#NEED CONVERSION OF PRACTICAL SALINITY TO ABSOLUTE SALINITY BEFORE USING THE FOLLOWING FUNCTION
#NEED CONVERSION OF IN-SITU TEMPERATURE TO CONSERVATIVE TEMPERATURE
psal <- gsw_SA_from_SP(finaldf$PSAL,finaldf$depth,finaldf$lon,finaldf$lat)
temp <- gsw_CT_from_t(psal,finaldf$TEMP,finaldf$depth)
rho_anomaly <- gsw_sigma0(psal,temp)

rhodf <- data.frame(Depth= finaldf$depth,
           rho_anomaly = rho_anomaly,  
           CHLA = finaldf$CHLA,
           TEMP = finaldf$TEMP,
           PSAL = psal,
           juld            = finaldf$juld,
           day             = month.day.year(finaldf$juld,c(1,1,1950))$day,
           month           = month.day.year(finaldf$juld,c(1,1,1950))$month,
           year            = month.day.year(finaldf$juld,c(1,1,1950))$year,
           lon             = finaldf$lon,
           lat             = finaldf$lat)
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
  
  #Extraction de l'anomalie de densité potentielle
  densitydf <- ExtractDensity(chladf, FloatInfo)
    
  subchladf <- subset(densitydf,select=c("rho_anomaly","juld","CHLA","lon","lat"))
    
  chladf <- ddply(subchladf,~juld,summarize,
                    rho_anomalymax = rho_anomaly[which.max(CHLA)],
                    maxvalue = CHLA[which.max(CHLA)],
                    integration = sum(CHLA),
                    lon=mean(lon),
                    lat=mean(lat))
    
    id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
    
    nc_close(ncfile)
    
    #Construction du data frame final a lieu ici
    data.frame(rho_anomalymax  = chladf$rho_anomalymax,
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

########################################################################### PROFILES DF (TFE)

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
  
  #Extraction de l'anomalie de densité potentielle
  densitydf <- ExtractDensity(chladf, FloatInfo)
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
  data.frame(density       = densitydf$rho_anomaly,
             juld            = densitydf$juld,
             chla           = densitydf$CHLA,
             day             = month.day.year(densitydf$juld,c(1,1,1950))$day,
             month           = month.day.year(densitydf$juld,c(1,1,1950))$month,
             year            = month.day.year(densitydf$juld,c(1,1,1950))$year,
             DOY             = as.integer(strftime(as.Date(densitydf$juld,origin = '1950-01-01'), format ="%j")),#Day Of Year
             #lon             = chladf$lon,
             #lat             = chladf$lat,
             Platform        = as.numeric(unique(id)),
             type            = "Argo")
})

#pour supprimer les doublons qui apparaissent... (note que l'on supprime peut-être trop.. à check)
profiledf <- unique(profiledf)
#Assign a profile number according to the (unique due to two decimals) julian day 
profiledf <- transform(profiledf,id=as.numeric(factor(juld)))

#éventuellement -> introduit un peu d'incertitude mais pas tant que cela
density_duplicate <- which(duplicated(profiledf$density))
profiledf <- profiledf[-density_duplicate,]


## PAS SUPER CODE ##
#profiles that we keep
keepid <- read.table("dat_juld.txt", header = TRUE, sep = "", numerals = c("warn.loss"))

#Reorder
profiledf <- profiledf[order(profiledf$juld, na.last=TRUE, decreasing=FALSE),]
n <- profiledf

#Remove unecessary profiles from the "database" profiledf
#float precision............
#all.equal
tmp <- profiledf[profiledf$DOY==keepid$DOY[1] &
            profiledf$Platform == keepid$Platform[1] &
              profiledf$year == keepid$year[1],]

for (i in 2:length(keepid$DOY)){
tmp2 <- profiledf[profiledf$DOY==keepid$DOY[i] &
                 profiledf$Platform == keepid$Platform[i] &
                 profiledf$year == keepid$year[i],]
tmp <- rbind(tmp,tmp2)
}

tmp <- transform(tmp,id=as.numeric(factor(juld)))#pas parfait mais bon
profiledf <- tmp

## END PAS SUPER CODE

###### TEST FILTRATION DONNÉES

filterpoints <- 10
  
# remove depth duplicate
testdf <- ldply(as.list(1:length(unique(profiledf$id))), function(i){
  tmp <- profiledf[profiledf$id==i,]
  if (length(tmp$chla) > filterpoints){
  meanchla <- SMA(tmp$chla,filterpoints)#Simple Moving Average (remove sharp peaks and smooth the global curve)
  tmp$chla <- meanchla
  data.frame(tmp = tmp)
  }
})

colnames(testdf) <- c("density", "juld","chla","day","month","year","DOY","Platform","type","id")

profiledf <- testdf

#########

profiledf2<-ddply(profiledf,~month, transform, season=1*(month %in% c(12,1,2))+
                    2*(month %in% c(3,4,5 ))+
                    3*(month %in% c(6,7,8 ))+
                    4*(month %in% c(9,10,11 )))

#plots (cfr. Navarro et al. 2013, figure 3)

# ggplot(profiledf, aes(x=density, y=chla, color=month, group=juld)) + geom_line() +
#   facet_grid(~year) + xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   coord_flip() + xlim(c(15.5,13.5)) 
# 
# ggplot(profiledf, aes(x=density, y=chla, color=month, group=juld)) + geom_line()+ 
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   coord_flip() + xlim(c(15.5,13.5))
# 
# ggplot(profiledf, aes(x=density, y=chla, color=month, group=juld)) +
#   geom_line() + facet_grid(~year) +
#   xlab("Depth (m)") + ylab("Chlorophyll a (kg/m³)") +
#   xlim(15, 13) + ylim(0,4) +
#   coord_flip() + scale_color_gradientn(colours = rainbow(5))

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

#same as DepthDataExtraction
ggplot(perioddf, aes(x=density, y=chla, group=juld)) +
  geom_line() + facet_grid(~group) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) + coord_flip() + scale_x_reverse()

ggplot(profiledf, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_viridis() +  scale_x_reverse()

ggplot(profiledf, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_gradientn(colours = rainbow(5)) + scale_x_reverse()

#for continuous scale
ggplot(profiledf2, aes(x=density, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_viridis() + scale_x_reverse()

#discrete case
profiledf2$season[profiledf2$season == 1] <- "Winter"
profiledf2$season[profiledf2$season == 2] <- "Spring"
profiledf2$season[profiledf2$season == 3] <- "Summer"
profiledf2$season[profiledf2$season == 4] <- "Autumn"

ggplot(profiledf2, aes(x=density, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_x_reverse()

profiledf2$month[profiledf2$month == 1] <- "Jan"
profiledf2$month[profiledf2$month == 2] <- "Feb"
profiledf2$month[profiledf2$month == 3] <- "Mar"
profiledf2$month[profiledf2$month == 4] <- "Apr"
profiledf2$month[profiledf2$month == 5] <- "May"
profiledf2$month[profiledf2$month == 6] <- "Jun"
profiledf2$month[profiledf2$month == 7] <- "Jul"
profiledf2$month[profiledf2$month == 8] <- "Aug"
profiledf2$month[profiledf2$month == 9] <- "Sep"
profiledf2$month[profiledf2$month == 10] <- "Oct"
profiledf2$month[profiledf2$month == 11] <- "Nov"
profiledf2$month[profiledf2$month == 12] <- "Dec"

profiledf2$month = factor(profiledf2$month, levels=c('Jan','Feb','Mar','Apr','May','Jun',
                                                 'Jul','Aug','Sep','Oct','Nov','Dec'))

ggplot(profiledf2, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~year) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_x_reverse()

#2017

tmp <- profiledf2[profiledf2$year == 2017,]

ggplot(tmp, aes(x=density, y=chla, color=season, group=juld)) +
  geom_line() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_x_reverse()

ggplot(profiledf2, aes(x=density, y=chla, color=month, group=juld)) +
  geom_line() + facet_grid(~season) +
  xlab("Potential density anomaly (kg/m³)") + ylab("Chlorophyll a (kg/m³)") +
  ylim(0,4) +
  coord_flip() + scale_color_gradientn(colours = rainbow(5)) + scale_x_reverse()
  coord_flip() + scale_color_viridis() + scale_x_reverse()

########################################################################################



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

argodf5 <- subset(argodf2, select=c("lon","lat","juld","rho_anomalymax","DOY"))
write.table(argodf5, file="rhoInput_DOY.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = TRUE)

argodf2$season[which(argodf2$season == 1)]<-"Winter"
argodf2$season[which(argodf2$season == 2)]<-"Spring"
argodf2$season[which(argodf2$season == 3)]<-"Summer"
argodf2$season[which(argodf2$season == 4)]<-"Autumn"

# Le fait de le définir en "factor" va faire qu'il va respecter l'ordre imposé
argodf2$season<-factor(argodf2$season,levels = c("Winter","Spring","Summer","Autumn"))

ggplot(argodf2, aes(x = rho_anomalymax, fill=factor(year))) + geom_density(alpha=.5) +
  facet_grid(.~season,scales = "free") + coord_flip() + xlim(c(16,12))

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
---------------------------------------------------------------------------------
  # Extracting data from WOD ----------------------------------------------------
---------------------------------------------------------------------------------
#Measured by CTD and OSD... -> En attendant de résoudre le problème de la base
#de données EMODnet (807 datasets dispos..) -> OK on garde la WOD??
  
  
### ATTENTION : Utilisateur doit mettre ses fichiers ship-based dans un 
### sous-répertoire sous le nom /Merged/For_rho si on utilise le dossier
### pour le calcul de l'anomalie de densité potentielle
wd <- getwd()

#Creation of a list of files
#tmp <- list.files(path = paste0(wd,"/Merged"), pattern="*.nc")
tmp <- list.files(path = paste0(wd,"/Merged/For_rho/"), pattern="*.nc")
# tmp <- c("wod_011529449O.nc")
# tmp <- c("wod_011523761O.nc")

#fonction qui va éliminer certains fichiers qui ne répondent pas à des critères
#que nous avons établis sur certaines bases
filedf <- ldply(as.list(tmp), function(file){
  bin <- 1#1 if criteria are ok, 0 otherwise
  ncfile <- nc_open(paste0(wd,"/Merged/For_rho/",file), write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  chla <- ncvar_get(ncfile,"Chlorophyll")
  nsamples <- length(chla)
  pressure <- ncvar_get(ncfile,"Pressure")
  maxpressure <- pressure[nsamples]
  lat <- ncvar_get(ncfile, "lat")
  lon <- ncvar_get(ncfile,"lon")
  nc_close(ncfile)
  
  #criteria
  if (nsamples < 3){
    bin <- 0
  }
  else if (maxpressure < 60){#attention seasonal effect -> ai pris une sorte de borne inf
    bin <- 0
  }
  else if (length(chla[chla < 0]) != 0){
    bin <- 0
  }
  #Pour virer le plateau continental nord-ouest (et aussi un peu sur l'extrême ouest)
  else if (lat > 44.7 | lon < 28.4 | (lat > 43.7 & lon < 31)){
    bin <- 0
  }
  data.frame(filename = file, bin = bin)
})

subfiledf <- subset(filedf, bin == 1)
subfiledf <- subset(subfiledf, select = "filename")
filelist <- as.list(subfiledf$filename) 

shipdf <- ldply(as.list(filelist),function(file){
  
  ncfile <- nc_open(paste0(wd,"/Merged/For_rho/",file), write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  chla <- ncvar_get(ncfile,"Chlorophyll")
  depth <- ncvar_get(ncfile,"Pressure")
  lat <- ncvar_get(ncfile,"lat")
  long <- ncvar_get(ncfile,"lon")
  time <- ncvar_get(ncfile,"time")
  id <- ncvar_get(ncfile,"wod_unique_cast")
  temp <- ncvar_get(ncfile,"Temperature")#Temperature in situ should be OK
  psal <- ncvar_get(ncfile,"Salinity")#Quid de psal? Need conversion from practical salinity to absolute ?
  psal <- gsw_SA_from_SP(psal,depth,long,lat)
  rho_anomaly <- gsw_sigma0(psal,temp)
  #portion de code ajoutée pour éliminer le problème des NA dans rho_anomalymax (dus à des bugs de capteurs)
  # index <- which(!is.na(rho_anomaly), arr.ind = TRUE)
  # index <- index[which.max(chla)]
  # rho_anomalymax <- rho_anomaly[index]
  # nc_close(ncfile)
  # data.frame(Chlorophyll = chla, Depth = depth, rho_anomalymax = rho_anomalymax, day = month.day.year(time,c(1,1,1770))$day,
  #            month = month.day.year(time,c(1,1,1770))$month, year = month.day.year(time,c(1,1,1770))$year,
  #            Latitude = lat, Longitude = long, Platform = id)
  ### Plus besoin, on ne perd que 7 fichiers
  nc_close(ncfile)
  data.frame(Chlorophyll = chla, Depth = depth, rho_anomaly = rho_anomaly, day = month.day.year(time,c(1,1,1770))$day,
             month = month.day.year(time,c(1,1,1770))$month, year = month.day.year(time,c(1,1,1770))$year,
             Latitude = lat, Longitude = long, Platform = id)
  
})

#Final data frame
shipdff <- ddply(shipdf,~Platform,summarize, 
                 rho_anomalymax = rho_anomaly[which.max(Chlorophyll)],
                 #rho_anomalymax = rho_anomalymax[which.max(Chlorophyll)],
                 maxvalue = Chlorophyll[which.max(Chlorophyll)],
                 integratedvalue = sum(Chlorophyll), day = day[which.max(Chlorophyll)],
                 month = month[which.max(Chlorophyll)], year = year[which.max(Chlorophyll)],
                 lon = Longitude[which.max(Chlorophyll)], lat = Latitude[which.max(Chlorophyll)],
                 type = "ship-based")#ship-based car pas d'infos précises sur WOD (CTD ou OSD en l'occurrence)

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

#On a ici de la marge pour supprimer des entrées en 'post-traitement' (par exemple des intégrations de
#chlorophylle beaucoup trop importantes voir d'autres critères..)

# #Data distribution
# ggplot(datadff, aes(depthmax)) + geom_density()#note qu'il y a bcp de zéro car il semble y avoir pas mal de sampling en surface... (intérêt?)
# 
# #Spatial Coverage
# myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
# ggmap(myMap) + geom_point(aes(x=Longitude, y=Latitude), data = datadff, alpha = .2)


#Merging of ARGO and ship-based samplings  -------------------------------------------
merged <- merge(argodf2,shipdf2,all=T)
merged <- subset(merged, select = c("rho_anomalymax","maxvalue", "integratedvalue", 
                                    "day", "month", "year", "lon", "lat",
                                    "Platform", "type", "season"))

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
#save animation
#gganimate(test, "general.gif")

### Using density as Y axis (in place of pressure) -------

ggplot(shipdf2, aes(x = rho_anomalymax)) +
  geom_density(alpha=.5) +coord_flip()+ xlim(c(15,11))+
  facet_grid(~season,scales = "free")

ggplot(argodf2, aes(x = rho_anomalymax)) +
  geom_density(alpha=.5)+coord_flip()+xlim(c(17.5,11))+
  facet_grid(~season,scales = "free")

ggplot(shipdf2, aes(x = rho_anomalymax)) +
  geom_density(alpha=.5)+coord_flip() + xlim(c(15,11))

ggplot(argodf2, aes(x = rho_anomalymax)) +
  geom_density(alpha=.5)+coord_flip()+xlim(c(17.5,11))

ggplot(shipdf2, aes(x = rho_anomalymax, fill=factor(year), frame = year)) +
  geom_density(alpha=.5)+coord_flip()+ xlim(c(15,11))+
  facet_grid(.~season,scales = "free")

z <- ggplot(argodf2, aes(x = rho_anomalymax, fill=factor(year), frame = year)) +
  geom_density(alpha=.5)+coord_flip()+xlim(c(16,11))+
  facet_grid(.~season,scales = "free")

gganimate(z,interval = 4)
# gganimate(z, "density_anomaly.gif")

#seasonal density plot for argo AND ship data
ggplot(merged, aes(x = rho_anomalymax, color = type)) + geom_density(alpha=.5) + coord_flip() +
  facet_grid(.~season, scales = "free") + xlim(c(16,11)) 

ggplot(merged, aes(x = rho_anomalymax, color = type)) + geom_density(alpha=.5) + coord_flip()

#################################################################################################
#################################################################################################
#################################################################################################

#Add downwelling photosynthetic radiation --------
#ATTENTION seulement 3 bouées sont équipées de capteur de PAR
filename <- c("6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

#Fonction similaire à ExtractDensity mais adapté à la PAR 
ExtractDensityPAR <- function(pardf, FloatInfo){
  
  tempdf <- ExtractVar("TEMP_ADJUSTED",FloatInfo)
  if (all(is.na(tempdf)) == TRUE) {
    tempdf <- ExtractVar("TEMP",FloatInfo)
  }
  
  psaldf <- ExtractVar("PSAL_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    psaldf <- ExtractVar("PSAL",FloatInfo)
  }
  
  #Extract matching rows of a data frame
  tempdf <- match_df(tempdf,pardf, on = c("depth", "juld"))
  psaldf <- match_df(psaldf,pardf, on = c("depth", "juld"))
  
  #Sub data frames
  
  subtempdf <- subset(tempdf, select = c("value","depth","aprofile", "alevel","juld","lon","lat"))
  colnames(subtempdf)[which(colnames(subtempdf)=="value")]<-"TEMP"
  subpsaldf <- subset(psaldf, select = c("value","depth","aprofile","juld"))
  colnames(subpsaldf)[which(colnames(subpsaldf)=="value")]<-"PSAL"
  subpardf <- subset(pardf, select = c("value","depth","aprofile", "alevel","juld"))
  colnames(subpardf)[which(colnames(subpardf)=="value")]<-"PAR"
  colnames(subpardf)[which(colnames(subpardf)=="aprofile")]<-"profilePAR"
  colnames(subpardf)[which(colnames(subpardf)=="alevel")]<-"levelPAR"
  joindf <- join(subtempdf,subpsaldf, by = c("depth", "aprofile", "juld"))
  colnames(joindf)[which(colnames(joindf)=="aprofile")]<-"profileTEMP/PSAL"
  colnames(joindf)[which(colnames(joindf)=="alevel")]<-"levelTEMP/PSAL"
  subpardf <- match_df(subpardf,joindf, on = c("depth", "juld"))
  finaldf <- join(subpardf,joindf, by = c("depth","juld"))
  
  #use of gsw package
  #NEED CONVERSION OF PRACTICAL SALINITY TO ABSOLUTE SALINITY BEFORE USING THE FOLLOWING FUNCTION
  #TEMP is OK because conservative temperature is the same as in-situ temperature (I guess)
  psal <- gsw_SA_from_SP(finaldf$PSAL,finaldf$depth,finaldf$lon,finaldf$lat)
  rho_anomaly <- gsw_sigma0(psal,finaldf$TEMP)
  
  rhodf <- data.frame(Depth= finaldf$depth,
                      rho_anomaly = rho_anomaly,  
                      PAR = finaldf$PAR,
                      TEMP = finaldf$TEMP,
                      PSAL = psal,
                      juld            = finaldf$juld,
                      day             = month.day.year(finaldf$juld,c(1,1,1950))$day,
                      month           = month.day.year(finaldf$juld,c(1,1,1950))$month,
                      year            = month.day.year(finaldf$juld,c(1,1,1950))$year,
                      lon             = finaldf$lon,
                      lat             = finaldf$lat)
}


argopardf<-ldply(as.list(filename),function(file){
  
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
  pardf <- ExtractVar("DOWNWELLING_PAR_ADJUSTED",FloatInfo)
  if (all(is.na(psaldf)) == TRUE) {
    pardf <- ExtractVar("DOWNWELLING_PAR",FloatInfo)
  }
  
  #Extraction de l'anomalie de densité potentielle
  pardf <- ExtractDensityPAR(pardf, FloatInfo)
  
  subpardf <- subset(pardf,select=c("rho_anomaly","juld","PAR","lon","lat"))
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  nc_close(ncfile)
  
  #Construction du data frame final a lieu ici
  data.frame(rho_anomaly = subpardf$rho_anomaly,
             PAR        = subpardf$PAR,
             juld            = subpardf$juld,
             day             = month.day.year(subpardf$juld,c(1,1,1950))$day,
             month           = month.day.year(subpardf$juld,c(1,1,1950))$month,
             year            = month.day.year(subpardf$juld,c(1,1,1950))$year,
             lon             = subpardf$lon,
             lat             = subpardf$lat,
             Platform        = unique(id),
             type            = "Argo") 
  
})

#################################################################################################
#################################################################################################
#################################################################################################

#Fluorescence data

#change directory
#Creation of a list of files
tmp <- list.files(path = paste0(wd,"/FLUO/"), pattern="*.nc")

#Vu la tête du contenu des fichiers on va attendre avant de faire le filtre.. 
#Va falloir trouver un moyen de ne pas virer le fichier s'il y a des valeurs négatives ou aberrantes
# -> fixer des seuils de validité et mettre des NA là où les données aberrantes se trouvent pour 
# jeter tout le fichier.. NOTE : Pourrait se faire également dans le filtre ici plus haut.

fluodf <- ldply(as.list(tmp), function(file){

  ncfile <- nc_open(paste0(wd,"/FLUO/",file), write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  
  lat <- ncvar_get(ncfile, "LATITUDE")#lat, lon et juld n'ont pas la même taille que les autres
  lon <- ncvar_get(ncfile,"LONGITUDE")
  juld <- ncvar_get(ncfile,"TIME")#days since 1950-01-01
  juld[juld < 3400] <- NA #faudra filtrer les fichiers pourris
  juld[juld > 25000] <- NA
  fluo <- ncvar_get(ncfile,"FLU2")
  #pres <- ncvar_get(ncfile,"PRES")
  #temp <- ncvar_get(ncfile,"TEMP")
  #psal <- ncvar_get(ncfile,"PSAL")
  #rep function? plutôt que boucle...
  # n <- length(fluo[,1])
  # m <- length(lat)
  # #creation of empty matrices for repeating data (lat,lon and juld)
  # nlat <- matrix(nrow = n, ncol = m)
  # nlon <- matrix(nrow = n, ncol = m)
  # njuld <- matrix(nrow = n, ncol = m)
  
#   for (i in 1:m){
#   nlat[,i] <- lat[i]
#   nlon[,i] <- lon[i]
#   njuld[,i] <- juld[i]
# }
  
  # psal <- gsw_SA_from_SP(psal,pres,lon,lat)
  # psal[psal > 35] <- NA #Black Sea : Salinité comprise entre X et X (get rid of outliers)
  # psal[psal < 5] <- NA
  nc_close(ncfile)
  
  data.frame(lat = lat, lon = lon, juld = juld,day = month.day.year(juld,c(1,1,1950))$day,
             month           = month.day.year(juld,c(1,1,1950))$month,
             year            = month.day.year(juld,c(1,1,1950))$year)
})

#Ajout des saisons                 
fluodf2<-ddply(fluodf,~month, transform, season=1*(month %in% c(12,1,2))+
                 2*(month %in% c(3,4,5 ))+
                 3*(month %in% c(6,7,8 ))+
                 4*(month %in% c(9,10,11 )))

fluodf2$season[which(fluodf2$season == 1)]<-"Winter"
fluodf2$season[which(fluodf2$season == 2)]<-"Spring"
fluodf2$season[which(fluodf2$season == 3)]<-"Summer"
fluodf2$season[which(fluodf2$season == 4)]<-"Autumn"

fluodf2$season<-factor(fluodf2$season,levels = c("Winter","Spring","Summer","Autumn"))

#Spatio-temporal coverage
myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
ggmap(myMap) + geom_point(aes(x=lon, y=lat, color = season), data = fluodf2, alpha = .4)  

#Seasonal Coverage
ggmap(myMap) + 
  geom_point(aes(x=lon, y=lat, color = season), data = fluodf2, alpha = .2) + 
  facet_grid(~season) 

#################################################################################################
#################################################################################################
#################################################################################################

#JUST FOR ARGO DATA FOR NOW
#Data export test on argodf2
argodf3 <- subset(argodf2, select=c("lon","lat","juld","rho_anomalymax","maxvalue","day","month","year"))
write.table(argodf3, file="rhoInput.txt", sep=" ", na = "NA", dec = ".", eol = "\r\n",
            row.names = FALSE, col.names = TRUE)

# write.table(argodf3, file="rhoInput.csv", sep=",", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

# write.table(argodf3, file="rhoInput.txt", sep="\t", na = "NA", dec = ".", eol = "\r\n",
#             row.names = FALSE, col.names = TRUE)

#TO DO: REORDER THE DATAFRAME FROM MIN TO MAX AND/OR return a table with min and max values for 
#each column

system("bash inputfile.sh")#execution of a shell script to have a 'clean' input file


