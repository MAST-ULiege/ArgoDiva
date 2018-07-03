
#################################################################
#   AUTHORS : BASED ON MARIN CORNEC AND LOUIS TERRATS SCRIPTS   #
#################################################################

#################################################################
#                             REMARKS                           #
# This script has been adapted and modified from the original   #
# one. It is not perfect (see the name of the downloaded        #
# NetCDF files, e.g. PROFILESXXXX.nc instead of XXXXXX.nc)      #
# Moreover, it may bug if your machine is not powerful enough,  #
# better not do it on a laptop                                  #
#################################################################


library(tidyverse)
library(e1071)
library(pbapply)
library(ncdf4)
library(ggplot2)
library(plyr)
library(reshape2)
library(chron)
library(viridis)
library(TTR)
library(minpack.lm)
library(nls2)
library(nlstools)
library(data.table)
library(gsw)
library(gridExtra)
library(car)
library(oce)
library(ggmap)
library(ggalt)
library(stringr)

# Area of study
lat_min <- 40
lat_max <- 46
lon_min <- 25
lon_max <- 43

# Temporal window
maximal_date <- 20180501e+06 # To exclude new float data which appear each day

# Paths => TO BE ADAPTED TO EACH MACHINE
path_to_data <- "/home/flo/R/TFE/tfe/DATA/ARGO/DATA" # Path to directory containing Argo data
path_to_profiles <- "/home/flo/R/TFE/tfe/DATA/ARGO/DATA/PROFILES"# Path to directory containing Argo profile data 

# Download the txt file (HUGE FILE BE CAREFUL IT TAKES TIME)
download.file('ftp://ftp.ifremer.fr/ifremer/argo/ar_index_global_prof.txt',
              "/home/flo/R/TFE/tfe/DATA/ARGO/DATA/ar_index_global_prof.txt",
              quiet = FALSE, mode = "w", cacheOK = TRUE)

# Import txt file (downloaded at ftp://ftp.ifremer.fr/ifremer/argo/argo_bio-profile_index.txt)
data_txt <- read.csv(file.path(path_to_data, "ar_index_global_prof.txt"), comment.char="#")

# Definition of temporal window (2012/01/02 : launch of first VIIRS mission)
data_txt <- data_txt[!(data_txt$date > maximal_date),]

# Create a column containing WMO name of the float
data_txt$wmo <- str_split(data_txt$file, "/", simplify = TRUE)[,2]

# Create a column containing profile name 
data_txt$profile <- substr(str_split(data_txt$file, "/", simplify = TRUE)[,4],1,nchar(str_split(data_txt$file, "/", simplify = TRUE)[,4])-3)

data <- data_txt[!(data_txt$latitude > lat_max | data_txt$latitude < lat_min),]
data <- data[!(data$longitude > lon_max | data$longitude < lon_min),]

#CLEAN NA'S ================> BE CAREFUL WITH THIS IF YOU WANT TO KEEP LINES WITH 'SOME' NA'S BUT NOT ALL !
data <- data[complete.cases(data),]

# Creation of a vector which contain all WMO codes, dac, name of the profile of floats of interest and the parameters
float_list <- setNames(data.frame(matrix(ncol = 3, nrow = nrow(data))), c("idd", "WMO", "profile"))
float_list$idd <- data$file
float_list$WMO <- data$wmo
float_list$YEAR <- substr(data$date, 1, 4)
float_list$profile <- data$profile
float_list$LAT <- data$latitude
float_list$LON <- data$longitude
float_list$TIME <- as.Date(substr(data$date, 1, 8), format = "%Y%m%d")

# Save of the float list in the txt format
write.table(float_list, file.path(path_to_data, "float_list.txt"), sep="\t")

already_downloaded <- grep(list.files(path_to_profiles, pattern = ".nc$", recursive = TRUE), pattern = paste(unique(float_list$WMO), collapse = "|"),
                           value = T)


# Variables to display progress of the downloading
tot <- dim(float_list)[1]
count <- length(already_downloaded)

for (j in 1:10) { # Downloading loop is launch 10 times (to download data that are not updated because of an error message)
  
  # Update of the list containing already downloaded profiles (update at each iterations)
  already_downloaded <- grep(list.files(path_to_profiles, pattern = ".nc$", recursive = TRUE), pattern = paste(unique(float_list$WMO), collapse = "|"), 
                             value = T)
  
  # Update of the variable to display the true progress
  count <- length(already_downloaded) 
  
  if (length(already_downloaded) != dim(float_list)[1]) { # If it misses profiles, relaunch the loop
    
    float_list_bis <- split(float_list, float_list$WMO)
    
    pblapply(float_list_bis, function(x) {
      
      sapply(x$idd, function(y, count) {
        # for (i in 1:dim(float_list)[1]) { # Downloading loop
        
        # If profile doesn't exist download it, otherwise nothing is done
        if (y %in% already_downloaded == FALSE) { 
          
          setwd(path_to_profiles)
          
          # Download of data and metadata
          try({ # try function to automatically keep up the for loop after an error message
            download.file(paste('ftp://ftp.ifremer.fr/ifremer/argo/dac/',y,sep = ""),
                          paste(path_to_profiles, str_split(y,"/", simplify = TRUE)[[4]], sep = ""), mode ="wb", quiet = TRUE)
            
          },
          
          silent = TRUE)
          
          # # Update of the variable to display the true progress
          # count <- count + 1
          # print(paste('ftp://ftp.ifremer.fr/ifremer/argo/dac/',y, "   n° = ", count, "/", tot, sep = "")) # Text showing the progress of the downloading
          
        }
        return(count)}, count)
    })
    
    if ((nrow(float_list) == length(already_downloaded)) & (all(already_downloaded == float_list$idd) == TRUE)) { 
      stop("downloading is completed (this error message occurs to stop the for loop, not to report an error)")
    } else {print("downloading in progress")}
  }
}
