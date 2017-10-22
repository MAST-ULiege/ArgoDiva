---
  title: "Data extraction from NetCDF files"
author: "Florian Ricour"
date: "October 19, 2017"
remark : "Only works for ARGO buoys with CHLA, CDOM and DOXY sensors with their respective name"
---
  
  # SECTION ONE : EXTRACTING BASIC DATA ----------------------------------------------  

library(ncdf4)
library(matrixStats)#Fast statistical analysis for matrices
#More complex libraries for plotting data
# library(ggplot2)
# library(ggvis)
# library(lattice)

#Setting the working directory
setwd("~/R/TFE/tfe")
filename <- "6901865_Mprof.nc"

#Opening the file in a open-only mode
ncfile <- nc_open(filename, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)

#Dimensions
N_PROF <- ncol(ncvar_get(ncfile,"PRES"))
N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))

#Data extraction
#PRESSURE
pres <- ncvar_get(ncfile,"PRES")
pres_adjusted <- ncvar_get(ncfile,"PRES_ADJUSTED")

#CHLOROPHYLL-A
chla <- ncvar_get(ncfile,"CHLA")
chla_adjusted <- ncvar_get(ncfile,"CHLA_ADJUSTED")

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
wd <- paste(wd,"/",nwd, sep ="")
dir.create(wd)
setwd(wd)

#NOT ajusted QC go from 0 to 4

#CHLOROPHYLL-A
for (i in 0:4){#There are 5 'RAW' QC
  #final temp
  if (i == 0){
    tmp <- (nchla_qc == i & chla) * chla
    tmp[tmp == 0] <- NA
    f <- list(tmp)
  }
  else{
    tmp <- (nchla_qc == i & chla) * chla
    tmp[tmp == 0] <- NA
    tmpi <- list(tmp)
    f <- c(f,tmpi)
  }
}

#Final plot
#Save plot in the new working directory
#Format png pour commencer
png(filename = paste(nwd,"_CHLA",sep =""))
plot(y = -pres, x = f[[4]], main = "Raw Values", xlab = "Chlorophyll-a (mg/m3)", ylab = "Pressure (decibar)", col = "red", cex = 0.25)#le point du premier 'plotting' doit avoir des valeurs sinon ça bugge
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[2]], col = "green", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
dev.off()

#Temperature
for (i in 0:4){
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

png(filename = paste(nwd,"_TEMP",sep =""))
plot(y = -pres, x = f[[2]], main = "Raw Values", xlab = "Temperature (°C)", ylab = "Pressure (decibar)", col = "green", cex = 0.25)
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
dev.off()

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

png(filename = paste(nwd,"_PSAL",sep =""))
plot(y = -pres, x = f[[2]], main = "Raw Values", xlab = "Salinity", ylab = "Pressure (decibar)", col = "green", cex = 0.25)
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
dev.off()

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

png(filename = paste(nwd,"_CDOM",sep =""))
plot(y = -pres, x = f[[1]], main = "Raw Values", xlab = "CDOM (ppb)", ylab = "Pressure (decibar)", col = "blue", cex = 0.25)
points(y = -pres, x = f[[2]], col = "green", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
dev.off()

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

#par(mar = c(5.1, 4.1, 4.1, 8.1), xpd =TRUE) #peut-être trouver un moyen d'ajuster la postion de la légende automatiquement
png(filename = paste(nwd,"_DOXY",sep =""))
plot(y = -pres, x = f[[2]], main = "Raw Values", xlab = "Dissolved Oxygen (micromole/kg)", ylab = "Pressure (decibar)", col = "green", cex = 0.25)
points(y = -pres, x = f[[1]], col = "blue", cex = 0.25)
points(y = -pres, x = f[[3]], col = "orange", cex = 0.25)
points(y = -pres, x = f[[4]], col = "red", cex = 0.25)
points(y = -pres, x = f[[5]], col = "black", cex = 0.25)
legend("bottomright",legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
#legend("topright", inset = c(-0.3,0), legend=c("QC = 0","QC = 1","QC = 2","QC = 3", "QC = 4"), pch = c(1,1), col = c("blue","green","orange","red","black"))
dev.off()

# SECTION FOUR : EXTRACTION OF CHLOROPHYLL-A FEATURES -----------------------------------------------------

# #valuemax <- apply(chla,2,max, na.rm = TRUE)
# valuemax[valuemax == -Inf] <- 0
# #integratedvalue <- apply(chla,2,sum,na.rm = TRUE)

#Clean chla matrix (all NA columns removed)
index <- which(is.na(chla[1,]) == FALSE)
for (i in 1:length(index)){
  if(i == 1){
    nchla <- chla[,index[i]]
  }
  else{
    tmp <- chla[,index[i]]
    nchla <- cbind(nchla,tmp)
  }
}
colnames(nchla) <- c(1:ncol(nchla)) 
 
valuemax <- colMaxs(nchla, na.rm = TRUE)
totvalue <- colSums2(nchla, na.rm = TRUE)
maxchla <- seq(0,0,length.out = ncol(nchla))
maxdepth <- seq (0,0,length.out = ncol(nchla))
for (i in 1 : ncol(nchla)){
  maxchla[i] <- which.max(nchla[,i])
  maxdepth[i] <- pres[maxchla[i],index[i]]
}

datachla <- data.frame(valuemax,totvalue,maxdepth)
colnames(datachla) <- c("Max CHLA","Total CHLA","Max Depth")
