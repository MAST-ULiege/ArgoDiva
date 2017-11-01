---
title: "Data extraction from NetCDF files"
author: "Florian Ricour"
date: "October 18, 2017"
remark : "Only works for ARGO buoys with CHLA, CDOM and DOXY sensors with their respective name"
---

# SECTION ONE : EXTRACTING BASIC DATA ----------------------------------------------  
  
library(ncdf4)
#More complex libraries for plotting data
# library(ggplot2)
# library(ggvis)
# library(lattice)

#Setting the working directory
setwd("~/R/TFE/tfe")
filename <- "6901866_Mprof.nc"

#Opening the file in a open-only mode
ncfile <- nc_open(filename, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)

#print(ncfile)#equivalent de ncdisp dans Matlab

#Dimensions
#N_PROF <- ncol(ncvar_get(ncfile,"STATION_PARAMETERS"))
N_PROF <- ncol(ncvar_get(ncfile,"PRES"))
N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))

#Data extraction
#PRESSURE
pres <- ncvar_get(ncfile,"PRES")
pres_adjusted <- ncvar_get(ncfile,"PRES_ADJUSTED")
#pres_FillValue <- ncatt_get(ncfile,"PRES_ADJUSTED_ERROR",'_FillValue')#Pris en compte dans le suppress_dimvals

#CHLOROPHYLL-A
chla <- ncvar_get(ncfile,"CHLA")
chla_adjusted <- ncvar_get(ncfile,"CHLA_ADJUSTED")
#chla_FillValue <- ncatt_get(ncfile,"CHLA_ADJUSTED_ERROR",'_FillValue')

#TEMPERATURE
temp <- ncvar_get(ncfile,"TEMP")
temp_adjusted <- ncvar_get(ncfile,"TEMP_ADJUSTED")

#SALINITY
psal <- ncvar_get(ncfile,"PSAL")
psal_adjusted <- ncvar_get(ncfile,"PSAL_ADJUSTED")

#CDOM
cdom <- ncvar_get(ncfile,"CDOM")
cdom_adjusted <- ncvar_get(ncfile,"CDOM_ADJUSTED")

#DISSOLVED OXYGEN
doxy <- ncvar_get(ncfile,"DOXY")
doxy_adjusted <- ncvar_get(ncfile,"DOXY_ADJUSTED")

# #Chlorophyll-A plot
# plot (x = chla, y = -pres, main = "Raw Values", xlab = "Chlorophyll-A (mg/m3)", ylab = "Pressure (decibar)", cex = 0.25)
# 
# #Temperature plot
# plot (x = temp, y = -pres, main = "Raw Values", xlab = "Temperature (°C)", ylab = "Pressure(decibar)", cex = 0.25)
# 
# #Salinity plot
# plot (x = psal, y = -pres, main = "Raw Values", xlab = "Salinity", ylab = "Pressure(decibar)", cex = 0.25)
# 
# #CDOM plot
# plot (x = cdom, y = -pres, main = "Raw Values", xlab = "CDOM (ppb)", ylab = "Pressure(decibar)", cex = 0.25)
# 
# #DOXY plot
# plot (x = doxy, y = -pres, main = "Raw Values", xlab = "Dissolved oxygen (micromole/kg)", ylab = "Pressure(decibar)", cex = 0.25)

# SECTION TWO : EXTRACTING QC -----------------------------------------------------  

pres_qc <- ncvar_get(ncfile,"PRES_QC")
pres_adjusted_qc <- ncvar_get(ncfile,"PRES_ADJUSTED_QC")

#CHLOROPHYLL-A
chla_qc <- ncvar_get(ncfile,"CHLA_QC")
chla_adjusted_qc <- ncvar_get(ncfile,"CHLA_ADJUSTED_QC")

#Split of character vectors in QC files
#Initialisation
nchla_qc <- as.numeric(unlist(strsplit(chla_qc[1], split="")))#nchla = chla_qc "clean" (-> n = new)

for (i in 2:N_PROF){#Il devrait y avoir moyen d'utiliser la fonction sapply pour faire la boucle...
  
  tmp <- as.numeric(unlist(strsplit(chla_qc[i], split="")))
  nchla_qc <- cbind(nchla_qc,tmp)
}

colnames(nchla_qc) <- c(1:N_PROF)

#Temperature
temp_qc <- ncvar_get(ncfile,"TEMP_QC")
temp_adjusted_qc <- ncvar_get(ncfile,"TEMP_ADJUSTED_QC")

#Split of character vectors in QC files
#Initialisation
ntemp_qc <- as.numeric(unlist(strsplit(temp_qc[1], split="")))

for (i in 2:N_PROF){
  
  tmp <- as.numeric(unlist(strsplit(temp_qc[i], split="")))
  ntemp_qc <- cbind(ntemp_qc,tmp)
}

colnames(ntemp_qc) <- c(1:N_PROF)

#Salinity
psal_qc <- ncvar_get(ncfile,"PSAL_QC")
psal_adjusted_qc <- ncvar_get(ncfile,"PSAL_ADJUSTED_QC")

#Split of character vectors in QC files
#Initialisation
npsal_qc <- as.numeric(unlist(strsplit(psal_qc[1], split="")))

for (i in 2:N_PROF){
  
  tmp <- as.numeric(unlist(strsplit(psal_qc[i], split="")))
  npsal_qc <- cbind(npsal_qc,tmp)
}

colnames(npsal_qc) <- c(1:N_PROF)

#CDOM
cdom_qc <- ncvar_get(ncfile,"CDOM_QC")
cdom_adjusted_qc <- ncvar_get(ncfile,"CDOM_ADJUSTED_QC")

#Split of character vectors in QC files
#Initialisation
ncdom_qc <- as.numeric(unlist(strsplit(cdom_qc[1], split="")))

for (i in 2:N_PROF){
  
  tmp <- as.numeric(unlist(strsplit(cdom_qc[i], split="")))
  ncdom_qc <- cbind(ncdom_qc,tmp)
}

colnames(ncdom_qc) <- c(1:N_PROF)

#DOXY
doxy_qc <- ncvar_get(ncfile,"DOXY_QC")
doxy_adjusted_qc <- ncvar_get(ncfile,"DOXY_ADJUSTED_QC")

#Split of character vectors in QC files
#Initialisation
ndoxy_qc <- as.numeric(unlist(strsplit(doxy_qc[1], split="")))

for (i in 2:N_PROF){
  
  tmp <- as.numeric(unlist(strsplit(doxy_qc[i], split="")))
  ndoxy_qc <- cbind(ndoxy_qc,tmp)
}

colnames(ndoxy_qc) <- c(1:N_PROF)

# SECTION THREE : PLOTS OF DATA ACCORDING QC VALUES -------------------------------

#set a new working directory to get plots according to the given filename
wd <- getwd()
nwd <- strsplit(filename,".nc")
nwd <- paste(wd,"/",nwd, sep ="")
dir.create(nwd)
setwd(nwd)

#NOT ajusted QC go from 0 to 4

#CHLOROPHYLL-A
nchla0 <- (nchla_qc == 0 & chla) * chla#chla data at QC = 0
nchla0[nchla0 == 0] <- NA 

nchla1 <- (nchla_qc == 1 & chla) * chla#on devrait adapter ça à une boucle
nchla1[nchla1 == 0] <- NA 

nchla2 <- (nchla_qc == 2 & chla) * chla
nchla2[nchla2 == 0] <- NA 

nchla3 <- (nchla_qc == 3 & chla) * chla
nchla3[nchla3 == 0] <- NA 

nchla4 <- (nchla_qc == 4 & chla) * chla
nchla4[nchla4 == 0] <- NA 

plot(y = -pres, x = nchla3, main = "Raw Values", xlab = "Chlorophyll-A (mg/m3)", ylab = "Pressure (decibar)", col = "red", cex = 0.25)
points(y = -pres, x = nchla4, main = "Raw Values", xlab = "Chlorophyll-A (mg/m3)", ylab = "Pressure (decibar)", col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 3", "QC = 4"), pch = c(1,1), col = c("red","black"))

#Temperature
ntemp0 <- (ntemp_qc == 0 & temp) * temp
ntemp0[ntemp0 == 0] <- NA 

ntemp1 <- (ntemp_qc == 1 & temp) * temp
ntemp1[ntemp1 == 0] <- NA 

ntemp2 <- (ntemp_qc == 2 & temp) * temp
ntemp2[ntemp2 == 0] <- NA 

ntemp3 <- (ntemp_qc == 3 & temp) * temp
ntemp3[ntemp3 == 0] <- NA 

ntemp4 <- (ntemp_qc == 4 & temp) * temp
ntemp4[ntemp4 == 0] <- NA 

plot(y = -pres, x = ntemp1, main = "Raw Values", xlab = "Temperature (°C)", ylab = "Pressure (decibar)", col = "green", cex = 0.25)
points(y = -pres, x = ntemp0, col = "blue", cex = 0.25)
points(y = -pres, x = ntemp2, col = "orange", cex = 0.25)
points(y = -pres, x = ntemp3, col = "red", cex = 0.25)
points(y = -pres, x = ntemp4, col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red", "black"))

#SECTION FOUR : TEST TO MAKE A ROUTINE OF WHAT HAVE BEEN DONE ABOVE --------------

#Temperature
for (i in 0:4){#There are 5 'RAW' QC
  #final temp
  if (i == 0){
    tmp <- (ntemp_qc == i & temp) * temp
    tmp[tmp == 0] <- NA
    f <- list(tmp)
  }
  else{
    tmp <- (ntemp_qc == i & temp) * temp
    tmp[tmp == 0] <- NA
    tmpi <- list(tmp)
    f <- c(f,tmpi)
  }
}

#Final plot
plot(y = -pres, x = f[[2]], main = "Raw Values", xlab = "Temperature (°C)", ylab = "Pressure (decibar)", col = "green", cex = 0.25)#le point du premier 'plotting' doit avoir des valeurs sinon ça bugge
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))

#Salinity
for (i in 0:4){
  if (i == 0){
    tmp <- (npsal_qc == i & psal) * psal
    tmp[tmp == 0] <- NA
    f <- list(tmp)
  }
  else{
    tmp <- (npsal_qc == i & psal) * psal
    tmp[tmp == 0] <- NA
    tmpi <- list(tmp)
    f <- c(f,tmpi)
  }
}

#Final plot
plot(y = -pres, x = f[[2]], main = "Raw Values", xlab = "Salinity", ylab = "Pressure (decibar)", col = "green", cex = 0.25)
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))

#CDOM
for (i in 0:4){
  if (i == 0){
    tmp <- (ncdom_qc == i & cdom) * cdom
    tmp[tmp == 0] <- NA
    f <- list(tmp)
  }
  else{
    tmp <- (ncdom_qc == i & cdom) * cdom
    tmp[tmp == 0] <- NA
    tmpi <- list(tmp)
    f <- c(f,tmpi)
  }
}

#Final plot
plot(y = -pres, x = f[[1]], main = "Raw Values", xlab = "CDOM (ppb)", ylab = "Pressure (decibar)", col = "blue", cex = 0.25)
points(y = -pres, x = f[[2]], col = "green", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))

#DOXY
for (i in 0:4){
  if (i == 0){
    tmp <- (ndoxy_qc == i & doxy) * doxy
    tmp[tmp == 0] <- NA
    f <- list(tmp)
  }
  else{
    tmp <- (ndoxy_qc == i & doxy) * doxy
    tmp[tmp == 0] <- NA
    tmpi <- list(tmp)
    f <- c(f,tmpi)
  }
}

#Final plot
plot(y = -pres, x = f[[2]], main = "Raw Values", xlab = "Dissolved Oxygen (micromole/kg)", ylab = "Pressure (decibar)", col = "green", cex = 0.25)
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
dev.off()
# #Save image
# wd <- getwd()#get working directory
# file <-"./img/doxy.png"
# if (file.exists(file)) stop(file, "already exists")
# dir.create(dirname(file), showWarnings = FALSE)
# png(file)
# dev.off()
