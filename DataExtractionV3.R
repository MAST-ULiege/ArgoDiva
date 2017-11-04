# SECTION ONE : EXTRACTING BASIC DATA FOR BIO-ARGO PROFILERS ------------------

library(ncdf4)
#More complex libraries for plotting data
library(ggplot2)
library(plyr)
#for melt function
library(reshape2)
#for julian date manipulation
library(chron)


#Seuls fichiers qui possèdent des données de CHLORO. Seuls deux d'entrent eux possèdent aussi du CDOM
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")
filename <- c("6901866_Mprof.nc")

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
  # On peut essayer de trouver méthode + simple pour chopper les indexs
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
  
  return(d=d) # Obligation de mettre un return? -> Pas sur, essaie...
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
  
  #SECTION THREE : USAGE OF PREVIOUS FUNCTIONS AND FINALIZATION OF OUR ARGO DATA FRAME ---------------
  
  CHLAL <- ExtractVar("CHLA",FloatInfo)
  
  CHLAL_ADJUSTED <- ExtractVar("CHLA_ADJUSTED",FloatInfo)
  
  #Subset of data frames
  subCHLA<-subset(CHLAL,select=c("depth","juld","value","qc","lon","lat"))
  colnames(subCHLA)[which(colnames(subCHLA)=="value")]<-"CHLA"
  
  subCHLA_ADJUSTED<-subset(CHLAL_ADJUSTED,select=c("depth","juld","value","qc","lon","lat"))
  colnames(subCHLA_ADJUSTED)[which(colnames(subCHLA_ADJUSTED)=="value")]<-"CHLA_ADJUSTED"
  
  CHLAF <- ddply(subCHLA,~juld,summarize,  qc          = qc[which.max(CHLA)], 
                                           depthmax    = depth[which.max(CHLA)],
                                           maxvalue    = CHLA[which.max(CHLA)],
                                           integration = sum(CHLA))
  
  CHLAF_ADJUSTED <- ddply(subCHLA_ADJUSTED,~juld,summarize,
                          qc = qc[which.max(CHLA_ADJUSTED)], 
                          depthmax = depth[which.max(CHLA_ADJUSTED)],
                          maxvalue = CHLA_ADJUSTED[which.max(CHLA_ADJUSTED)],
                          integration = sum(CHLA_ADJUSTED), 
                          lon=mean(lon),
                          lat=mean(lat))
  
  #RMQ : Comme on s'y attendait, la valeur de depthmax ne varie pas avec l'utilisation des données ajustées.
  #On garde les lignes de commande associées car elles auront leur utilité plus tard (max value et intégration)
  
  #Comme on a des données ajustées pour la chlorophylle pour les 4 fichiers, je n'utilise plus que les 
  #'ADJUSTED' values mais si ça n'était pas le cas, on aurait dû mettre une condition dès le début...

  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
  data.frame(depthmax        = CHLAF_ADJUSTED$depthmax,
             qc              = CHLAF_ADJUSTED$qc,
             maxvalue        = CHLAF_ADJUSTED$maxvalue,
             integratedvalue = CHLAF_ADJUSTED$integration,
             juld            = CHLAF_ADJUSTED$juld,
             day             = month.day.year(CHLAF_ADJUSTED$juld,c(1,1,1950))$day,
             month           = month.day.year(CHLAF_ADJUSTED$juld,c(1,1,1950))$month,
             year            = month.day.year(CHLAF_ADJUSTED$juld,c(1,1,1950))$year,
             lon             = CHLAF_ADJUSTED$lon,
             lat             = CHLAF_ADJUSTED$lat,
             Platform        = unique(id))
  
})

#SECTION FOUR : Three visualization of distribution density for different period (showing different approach)
ggplot(argodf, aes(x = depthmax, fill=factor(month))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))

# function "cut" provide a partioning based on variable month in 4 class
ggplot(argodf, aes(x = depthmax, fill=cut(month,4))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))

# --> Il semble qu'il y ait bien un cycle saisonnier marqué
# Maintenant je voudrais définir les saisons à mon gout et voir si ça correspond d'année en année

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

ggplot(argodf2, aes(x = depthmax, fill=factor(year))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(.~season,scales = "free")

# --> interessant, pour la bouée 6901866, la distribution en été et fort semblable pour 2015 et 2017 et bien différentes pour 2016...
# --> Dû à la variabilité spatiale ? 

ggplot(argodf2, aes(x = depthmax, fill=factor(year))) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(.~cut(lon,2),scales = "free")

#SECTION FIVE : SPATIAL COVERAGE OF ARGO PROFILERS ----------------------------------

### --> Ici j'ai plutot ajouté lon et lat dans le dataframe de base .. mais bon ça se discute .. 
#creation of a data frame for ARGO's trajectory
# trajdf <- ldply(as.list(filename), function(file){
#   ncfile <- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
#   id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
#   lat <- ncvar_get(ncfile, "LATITUDE")
#   long <- ncvar_get(ncfile, "LONGITUDE")
#   juld <- ncvar_get(ncfile,"JULD")
#   
#   nc_close(ncfile)
#   data.frame( Latitude   = lat,
#               Longitude  = long,
#               Platform   = id,
#               juld       = juld,
#               day        = month.day.year(juld,c(1,1,1950))$day,
#               month      = month.day.year(juld,c(1,1,1950))$month,
#               year       = month.day.year(juld,c(1,1,1950))$year
#   )
# })


#Trajectories
#ggplot(trajdf,aes(x=Longitude,y=Latitude,color=Platform)) + geom_path() 

#Map libraries
library(ggmap)
library(ggalt)

#Black Sea coordinates
bs <- c(26.5,40,43,46)
myMap <-get_map(location=bs, source="google",  crop=FALSE)
ggmap(myMap) +
  geom_point(aes(x=lon, y=lat, color=factor(season) ),
            data = argodf2, alpha = .8)+facet_grid(year~Platform)


---------------------------------------------------------------------------------
  #SECTION SIX : Data from WOD ----------------------------------------------------
---------------------------------------------------------------------------------
  #Measured by CTD and OSD... -> En attendant de résoudre le problème de la base
  #de données EMODnet (807 datasets dispos..)
  
wd <- getwd()
setwd(paste0(wd,"/Merged"))
#Creation of a list of files
tmp <- list.files(pattern="*.nc")
#tmp <- c("wod_011527263O.nc","wod_011529213O.nc")
datadf <- ldply(as.list(tmp),function(file){
  
  ncfile <- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  datachla <- ncvar_get(ncfile,"Chlorophyll")
  datadepth <- ncvar_get(ncfile,"Pressure")
  datalat <- ncvar_get(ncfile,"lat")
  datalong <- ncvar_get(ncfile,"lon")
  datadate <- ncvar_get(ncfile,"date")
  nc_close(ncfile)#sinon ça buggera aux alentours du 1000ème fichier..
  data.frame(Chlorophyll = datachla, Depth = datadepth, Latitude = datalat, Longitude = datalong, Date = datadate)
  
})

#Final data frame (dff)
datadff <- ddply(datadf,~Date,summarize, 
                 depthmax = Depth[which.max(Chlorophyll)],
                 maxvalue = Chlorophyll[which.max(Chlorophyll)],
                 integration = sum(Chlorophyll), Latitude = Latitude[which.max(Latitude)],
                 Longitude = Longitude[which.max(Longitude)])

#Data distribution
ggplot(datadff, aes(depthmax)) + geom_density()#note qu'il a bcp de zéro car il semble y avoir pas mal de sampling en surface... (intérêt?)

#Spatial Coverage
myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
ggmap(myMap) + geom_point(aes(x=Longitude, y=Latitude), data = datadff, alpha = .2)

#SECTION SEVEN : MERGING ARGO AND 'DISCRETE' SAMPLINGS -------------------------------------------
argodepthmax <- melt(argodf$depthmax)#one column data frame
discretedepthmax <- melt(datadff$depthmax)
merged <- rbind(argodepthmax, discretedepthmax)
colnames(merged)[which(colnames(merged)=="value")]<-"depthmax"

ggplot(merged, aes(depthmax)) + geom_density()+xlim(c(0,80))
