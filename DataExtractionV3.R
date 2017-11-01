# SECTION ONE : EXTRACTING BASIC DATA FOR BIO-ARGO PROFILERS ------------------

library(ncdf4)
#More complex libraries for plotting data
library(ggplot2)
library(plyr)
#for melt function
library(reshape2)

#Seuls fichiers qui possèdent des données de CHLORO. Seuls deux d'entrent eux possèdent aussi du CDOM
filename <- c("6900807_Mprof.nc","6901866_Mprof.nc","7900591_Mprof.nc","7900592_Mprof.nc")

#Construction of an ARGO-only dataframe
argodf<-ldply(as.list(filename),function(file){

#Opening the file in a open-only mode
ncfile <- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)

#Dimensions
N_PROF <- ncol(ncvar_get(ncfile,"PRES"))
N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))

juld <- ncvar_get(ncfile,"JULD")
pres <- ncvar_get(ncfile,"PRES")

#SECTION TWO : one function to extract variable and return a dataframe ---------------

ExtractVar<-function(Var){
  # This function should return a dataframe for variable with value, qc, iprofile and ilevel
  lvar             <- ncvar_get(ncfile,Var)
  lvar_qc          <- ncvar_get(ncfile,paste0(Var,"_QC"))
  lvar_qctab<-llply(lvar_qc,function(qcstring){
    as.numeric(unlist(strsplit(qcstring,split="")))
  })
  lvar_qctab<-do.call(cbind,lvar_qctab)
  
  # making dataframes, removing the NANs                 #On peut essayer de trouver méthode + simple pour chopper les indexes
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
  
  d$juld<-juld[d$aprofile]
  return(d=d)#Obligation de mettre un return?
}  

#SECTION THREE : USAGE OF PREVIOUS FUNCTIONS AND FINALIZATION OF OUR ARGO DATA FRAME ---------------

CHLAL <- ExtractVar("CHLA")
CHLAL_ADJUSTED <- ExtractVar("CHLA_ADJUSTED")

#Subset of data frames
subCHLA<-subset(CHLAL,select=c("depth","juld","value","qc"))
colnames(subCHLA)[which(colnames(subCHLA)=="value")]<-"CHLA"

subCHLA_ADJUSTED<-subset(CHLAL_ADJUSTED,select=c("depth","juld","value","qc"))
colnames(subCHLA_ADJUSTED)[which(colnames(subCHLA_ADJUSTED)=="value")]<-"CHLA_ADJUSTED"

CHLAF <- ddply(subCHLA,~juld,summarize, qc = qc[which.max(CHLA)], 
               depthofmaxCHLA = depth[which.max(CHLA)],
               maxvalue = CHLA[which.max(CHLA)],
               integration = sum(CHLA))

CHLAF_ADJUSTED <- ddply(subCHLA_ADJUSTED,~juld,summarize, qc = qc[which.max(CHLA_ADJUSTED)], 
                        depthofmaxCHLA = depth[which.max(CHLA_ADJUSTED)],
                        maxvalue = CHLA_ADJUSTED[which.max(CHLA_ADJUSTED)],
                        integration = sum(CHLA_ADJUSTED))

#RMQ : Comme on s'y attendait, la valeur de depthmax ne varie pas avec l'utilisation des données ajustées.
#On garde les lignes de commande associées car elles auront leur utilité plus tard (max value et intégration)

#Comme on a des données ajustées pour la chlorophylle pour les 4 fichiers, je n'utilise plus que les 
#'ADJUSTED' values mais si ça n'était pas le cas, on aurait dû mettre une condition dès le début...

#Construction du data frame final a lieu ici
data.frame(depthmax = CHLAF_ADJUSTED$depthofmaxCHLA,
           qc      = CHLAF_ADJUSTED$qc,
           maxvalue = CHLAF_ADJUSTED$maxvalue,
           integratedvalue = CHLAF_ADJUSTED$integration)

})

#SECTION FOUR : First attempt to get some kind of statistical distribution for ARGO data only ! ------------

ggplot(argodf, aes(x = depthmax)) + geom_bar() + xlim(c(2, 80)) #+ ylim(c(0,10))#limite des deux mètres (cfr. CTD)

ggplot(argodf, aes(x = depthmax, y = maxvalue)) + geom_point() + xlim(c(2, 80))

# ggplot(subset(argodf,qc==2),aes(x=depthmax))+ geom_bar()
# 
# ggplot(subset(argodf,qc==2),aes(x=depthmax, y = maxvalue))+ geom_point()
# 
# ggplot(subset(argodf,qc==3),aes(x=depthmax))+ geom_bar()
# 
# ggplot(subset(argodf,qc==3),aes(x=depthmax, y = maxvalue))+ geom_point()

plot(density(x = argodf$depthmax))

ggplot(argodf, aes(depthmax)) + geom_density()

#SECTION FIVE : SPATIAL COVERAGE OF ARGO PROFILERS ----------------------------------
#creation of a data frame for ARGO's trajectory
trajdf <- ldply(as.list(filename), function(file){
  
  ncfile <- nc_open(file, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)
  id <- ncvar_get(ncfile, "PLATFORM_NUMBER")
  lat <- ncvar_get(ncfile, "LATITUDE")
  long <- ncvar_get(ncfile, "LONGITUDE")
  data.frame(Latitude = lat, Longitude = long, Platform = id)
    
})

#Trajectories
#ggplot(trajdf,aes(x=Longitude,y=Latitude,color=Platform)) + geom_point() 

#Map libraries
library(ggmap)
library(ggalt)

#Black Sea coordinates
bs <- c(26.5,40,43,46)
myMap <-get_map(location=bs, source="google", maptype="satellite", crop=FALSE)
ggmap(myMap) + geom_point(aes(x=Longitude, y=Latitude, colour = Platform), data = trajdf, alpha = .2)



