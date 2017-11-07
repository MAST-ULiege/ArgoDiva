library(ncdf4)
#More complex libraries for plotting data
library(ggplot2)
library(plyr)
#for melt function
library(reshape2)
#for julian date manipulation
library(chron)

# SECTION ONE : EXTRACTING BASIC DATA FOR BIO-ARGO PROFILERS ------------------

#Seuls fichiers qui possèdent des données de CHLORO. Seuls deux d'entrent eux possèdent aussi du CDOM
#filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")
filename <- "6901866_Mprof.nc"

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
                          lon=mean(lon),
                          lat=mean(lat))#mean car la bouée a bougé sur la journée
  
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  
  #Construction du data frame final a lieu ici
  data.frame(depthmax        = chladf$depthmax,
             qc              = chladf$qc,
             maxvalue        = chladf$maxvalue,
             integratedvalue = chladf$integration,
             juld            = chladf$juld,#pour moi, on peut virer juld ici
             day             = month.day.year(chladf$juld,c(1,1,1950))$day,
             month           = month.day.year(chladf$juld,c(1,1,1950))$month,
             year            = month.day.year(chladf$juld,c(1,1,1950))$year,
             lon             = chladf$lon,
             lat             = chladf$lat,
             Platform        = unique(id),
             type            = "Argo")
  
})


#SECTION FOUR : Three visualization of distribution density for different period (showing different approach)
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
myMap <-get_map(location=bs, source="google", crop=FALSE)
ggmap(myMap) +
  geom_point(aes(x=lon, y=lat, color=factor(season) ),
            data = argodf2, alpha = .8)+facet_grid(year~Platform)

---------------------------------------------------------------------------------
  #SECTION SIX : Data from WOD ----------------------------------------------------
---------------------------------------------------------------------------------
  #Measured by CTD and OSD... -> En attendant de résoudre le problème de la base
  #de données EMODnet (807 datasets dispos..)
  
  
  
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
  nc_close(ncfile)
  
  #criteria
  if (nsamples < 3){
    bin <- 0
  }
  else if (maxpressure < 20){#test for 20m, attention seasonal effect...
    bin <- 0
  }
  else if (length(chla[chla < 0]) != 0){
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
# ggplot(datadff, aes(depthmax)) + geom_density()#note qu'il a bcp de zéro car il semble y avoir pas mal de sampling en surface... (intérêt?)
# 
# #Spatial Coverage
# myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
# ggmap(myMap) + geom_point(aes(x=Longitude, y=Latitude), data = datadff, alpha = .2)


#SECTION SEVEN : MERGING ARGO AND SHIP-BASED SAMPLINGS -------------------------------------------
merged <- merge(argodf2,shipdf2,all=T)
merged <- subset(merged, select = c("depthmax", "maxvalue", "integratedvalue", 
                                    "day", "month", "year", "lon", "lat",
                                    "Platform", "type", "season"))

ggplot(shipdf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(.~season,scales = "free")

ggplot(argodf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))+
  facet_grid(.~season,scales = "free")

ggplot(shipdf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))

ggplot(argodf2, aes(x = depthmax)) +
  geom_density(alpha=.5)+ xlim(c(2, 80))+coord_flip()+xlim(c(80,0))
