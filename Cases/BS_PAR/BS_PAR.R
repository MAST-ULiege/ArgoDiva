# Generation of the Black Sea Oxygen content files from Argo data for CMEMS OMI
# A. Capet 02/07/2018 - acapet@uliege.be

# Load package and function sources
source("../../Utils/ArgoLoad.R")
# Overwrites default ArgoSelect function
source("../../Utils/ArgoSelect_CMEMS.R")

datacase <- "CORIOLIS" #"CORIOLIS" #"CMEMS_HISTORY" #  

################################
# Criterium for Argo selection #
################################

source(paste0( "BSPARLoad_", datacase,".R" ) )

########################
# Function definitions #
########################
# To add new function, use :  
source("ArgoVertFunctions_BSPAR.R")
source("../../Utils/ArgoCompleteBIS.R")
source("../../Utils/ArgoExtractForCastedDF.R")

########################
# Complete (by level)  #
#      and / or        #
# Extract (by profile) #
########################

ddi  <- dcast(fulldf,qc+depth+juld+lon+lat+Platform+day+month+year~variable,fun.aggregate = mean, na.rm=T)

ddin <- subset(ddi,!is.na(PAR))
ddim <- melt(ddin,id.vars = c("qc","depth","juld","lon","lat","Platform","day","month","year"))

# No need for in situ functions here
#flistl <- list("InSituDens")
#fdf   <- ArgoCompleteBIS(ddim, flistl)
#ddi<-dcast(fdf,qc+depth+juld+lon+lat+Platform+day+month+year~variable,fun.aggregate = mean, na.rm=T)
#ddin<-subset(ddi,!is.na(DOXY) & !is.na(RHO))
#ddim<-melt(ddin,id.vars = c("qc","depth","juld","lon","lat","Platform","day","month","year"))

flist   <- list("PAR2bands")

finaldf <- ArgoExtractForCastDF(ddin, flist )
finaldf <- melt(finaldf, id.vars = c("qc","juld","lon","lat","Platform","day","month","year"))

ggplot(subset(finaldf, !is.na(value)), aes(x=value))+
  geom_density()+facet_wrap(~variable, scales = "free")

finaldffrop<-dcast(finaldf,qc+juld+lon+lat+Platform+day+month+year~variable, fun.aggregate = mean, na.rm=TRUE)

summary(finaldffrop)

#################

zforp <- seq(0.5,100,.5)

sumfinaldf <- ddply(finaldffrop, .(month),summarize, meanP0=mean(PAR0, na.rm=T))

PARCOMP<-ddply(sumfinaldf,.(month), function(dsub){
 data.frame(PAR=c( Att_2band(zforp,dsub$meanP0,0.23,0.48,0.137),
                   Att_2band(zforp,dsub$meanP0*0.5,0.23,0.48,0.137),
                   Att_2band(zforp,dsub$meanP0*1.5,0.23,0.48,0.137)),
                            depth=c(zforp,zforp,zforp), 
            group=c(rep(1,length(zforp)),rep(2,length(zforp)) ,rep(3,length(zforp))  ),
              month=dsub$month)})

ggplot(ddin,aes(x=PAR, y=-depth, color=factor(month)))+
  geom_point()+xlim(c(0,2000))+ylim(c(-100,0))+facet_wrap(~month, scales = "free_x")+
  geom_line(aes(group=group),data = PARCOMP, color="black")




