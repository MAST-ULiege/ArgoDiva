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
#TEMP is OK because conservative temperature is the same as in-situ temperature (I guess)
psal <- gsw_SA_from_SP(finaldf$PSAL,finaldf$depth,finaldf$lon,finaldf$lat)
rho_anomaly <- gsw_sigma0(psal,finaldf$TEMP)

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
               lon             = chladf$lon,
               lat             = chladf$lat,
               Platform        = unique(id),
               type            = "Argo") 
  
})


#Visualization of distribution density for different period (showing different approach)
ggplot(argodf, aes(x = rho_anomalymax)) + geom_density(alpha=.5) + coord_flip()

# Nasty code to get a user-defined season
argodf2<-ddply(argodf,~month, transform, season=1*(month %in% c(12,1,2))+
                 2*(month %in% c(3,4,5 ))+
                 3*(month %in% c(6,7,8 ))+
                 4*(month %in% c(9,10,11 )))

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
  else if (maxpressure < 20){#attention seasonal effect -> ai pris une sorte de borne inf
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

