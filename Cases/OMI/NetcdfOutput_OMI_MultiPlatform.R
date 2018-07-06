NetcdfOutput_OMI_MP <- function(df,varname){
  
  # -- Path here to save the netcdf file
  path <- paste0('./Output_',Sys.Date())
  
  # -- OMI file name : let's build the standard name of your future file
  area            <- 'bs'             # or med or nws, ...
  short_omi_name  <- 'oc'             # 
  first_occurence <- '20110101'       # yyyymmdd
  prod_time       <- paste0('P',gsub("-",'',as.character(Sys.Date())))      # P + yourdate yyyymmdd
  other_info      <- varname               # ohc case ...
  
  # -- TIME info
  TIMEdf <- data.frame(origin   = as.Date('1950-01-01 00:00:00'),
                       calendar = 'gregorian',
                       units    = 'days since 1950-01-01 00:00:00',
                       name     = 'time')

  # -- VAR infos
  VOCdf   <- data.frame(long_name     = "Oxygen inventory in the whole water column",
                       units         = 'mol/m^2',
                       name          = 'oxygen_content',
                       standard_name = 'oxygen_content',
                       FillValue     = -9999.0,
                       Comment       = 'Vertically integrated oxygen concetration',
                       minsc         = 0, 
                       maxsc         = 40,
                       ncname        = 'oxygen_inventory')
  rownames(VOCdf)<-"VOC"
  
  Z20df   <- data.frame(long_name     = "Oxygen Penetration Depth",
                       units         = 'm',
                       name          = 'oxycline_depth',
                       standard_name = 'oxycline_depth',
                       FillValue     = -9999.0,
                       Comment       = 'Depth at which [O2]<20µM', 
                       minsc         = 170, 
                       maxsc         = 0,
                       ncname        = 'oxygen_penetration_depth')
  rownames(Z20df)<-"Z20"
  
  R20df   <- data.frame(long_name     = "Oxygen Penetration Density",
                        units         = 'kg/m3',
                        name          = 'oxycline_density',
                        standard_name = 'oxycline_density',
                        FillValue     = -9999.0,
                        Comment       = 'Potential density anomaly at which [O2]<20µM',
                        minsc         = 16.5, 
                        maxsc         = 14.5,
                        ncname        = 'oxygen_penetration_density')
  
  rownames(R20df)<-"R20"
  
  allvardf <- rbind(VOCdf,Z20df,R20df)

  vardf <- allvardf[varname,]  
  
  ### ### ### ### #### 
  # Data Preparation #
  ### ### ### ### #### 

  df$juld <- round(df$juld)
  df      <- subset(df,select = c("juld","lon","lat","Platform",varname))
  colnames(df)[which(colnames(df)==varname)]<-"value"
  
  timevals <- seq(min(df$juld),max(df$juld))

  # -- Create platform list
  Platformlist <- unique(subset(df, !is.nan(value),select='Platform'))
  
  longdf <- data.frame(time = timevals, 
                       mean = NA,
                       stde  = NA)
  
  for (p in Platformlist$Platform) {
    print(p)
    platdf <- subset(df,Platform==p)
    longdf[match(platdf$juld ,longdf$time), p ] <- platdf$value
  }

  timespan <- 90 # days for a smoothed average and std error estimates  
  ##########################################
  # !!! Different option for stde definition !!! 
  # Pick your choice and comment the rest
  ##########################################
  # 1. simplest standard error on Mean estimate. -> problem due to sd estimation from very few floats (eg. 1)
  errorfromneighb <- function(neihgbdf) {
      sd( neihgbdf$value,na.rm = TRUE) / sqrt(sum(!is.na(neihgbdf$value) ))
  }
  
  # 2. Second option, the standard deviation used for standard error on mean estimate is estimated from the global datasets (not only the localo one)
  fullsd <- sd( df$value,na.rm = TRUE)
  errorfromneighb <- function(neihgbdf) {
  fullsd  / sqrt(sum(!is.na(neihgbdf$value) ))
  }
  
  for (d in longdf$time){
    if ( (d-min(longdf$time))<timespan/2 | (max(longdf$time)-d< timespan/2) ){longdf[which(longdf$time ==d), c('mean','stde')]<-NA}
    else{
    neihgbdf  <-  subset(df, abs(juld-d)<timespan/2)
    longdf[which(longdf$time ==d),'mean'] <- mean( neihgbdf$value,na.rm = TRUE)
    longdf[which(longdf$time ==d),'stde'] <- errorfromneighb(neihgbdf)
    }
  }

  ##########
  # Figure #
  ##########
  source("OMI_figure.R")
  figname <- paste0(area,'_omi_',short_omi_name,'_area_averaged_',vardf[,'ncname'],'_',first_occurence,'_',prod_time,'.pdf')
  OMI_figure(longdf, TIMEdf, vardf, figname)
  
  ###############
  # NC Creation #
  ###############
  ncfname <- paste0(area,'_omi_',short_omi_name,'_area_averaged_',vardf[,'ncname'],'_',first_occurence,'_',prod_time,'.nc')
  
  # -- Create the dimensions of the netcdf file
  timedim <- ncdim_def( as.character(TIMEdf$name)  ,
                        as.character(TIMEdf$units) ,
                        vals  = timevals           ,
                        unlim = TRUE)

  # -- Create vars
  count<-1
  defVar<-list()
  for (v in rownames(vardf)){
    for (p in Platformlist$Platform) {
    tmp<- ncvar_def(name     = paste0(vardf[v,'name'],'_',p),
              units    = as.character(vardf[v,'units']),
              dim      = list(timedim),
              missval  = vardf[v,'FillValue'], 
              longname = as.character(vardf[v,'long_name']))    
      
      defVar[[count]] <- tmp
      count<- count+1
    }
    
   tmp<- ncvar_def(name     = paste0(vardf[v,'name'],'_','mean'),
                      units    = as.character(vardf[v,'units']),
                      dim      = list(timedim),
                      missval  = vardf[v,'FillValue'], 
                      longname = as.character(vardf[v,'long_name']))    
   defVar[[count]] <- tmp
   count<- count+1
   
   tmp<- ncvar_def(name     = paste0(vardf[v,'name'],'_','std'),
                   units    = as.character(vardf[v,'units']),
                   dim      = list(timedim),
                   missval  = vardf[v,'FillValue'], 
                   longname = as.character(vardf[v,'long_name']))    
   defVar[[count]] <- tmp
   count<- count+1
  }
  

  # -- Create netcdf file  
  if (file.exists(ncfname)) file.remove(ncfname)
  ncout   <- nc_create(ncfname,defVar,force_v4=T)
  
  #############
  # ATRIBUTES #
  #############
  
  # -- Variable Attributes 
  #  - time
  ncatt_put(ncout,"time","standard_name","time")
  ncatt_put(ncout,"time","unit",TIMEdf$units)
  ncatt_put(ncout,"time","axis",'T')
  ncatt_put(ncout,"time","calendar",TIMEdf$calendar)

  for(v in rownames(vardf)){
    for (p in Platformlist$Platform) {  
      ncatt_put(ncout,paste0(vardf[v,'name'],'_',p), "standard_name", vardf[v,'standard_name'])
      ncatt_put(ncout,paste0(vardf[v,'name'],'_',p), "comment"      , paste0( vardf[v,'comment'], ' as measured from ARGO ', p))
      ncvar_put(ncout,paste0(vardf[v,'name'],'_',p),longdf[,p])
    }
    ncatt_put(ncout,paste0(vardf[v,'name'],'_','mean'), "comment"      , paste0( vardf[v,'comment'],' ',timespan, ' days window average'))
    ncvar_put(ncout,paste0(vardf[v,'name'],'_','mean'), longdf[,'mean'])
    ncatt_put(ncout,paste0(vardf[v,'name'],'_','std'), "comment"      , paste0( vardf[v,'comment'],' ',timespan, ' days window standard error'))
    ncvar_put(ncout,paste0(vardf[v,'name'],'_','std'), longdf[,'stde'])
  }
  

#  Global Attributes
  ncatt_put(ncout,0,'title'       ,'Black Sea Oxygen Inventory' )# to be consistent with PIT
  ncatt_put(ncout,0,'institution' ,'ULiege')
  ncatt_put(ncout,0,'references'  ,'http://marine.copernicus.eu') #please do not change
  ncatt_put(ncout,0,'Conventions' ,'CF-1.7')
  ncatt_put(ncout,0,'credit'      ,'E.U. Copernicus Marine Service Information (CMEMS)' )#please do not change
  ncatt_put(ncout,0,'contact'     ,'servicedesk.cmems@mercator-ocean.eu'       )         #please do not change
  ncatt_put(ncout,0,'source'      ,'INSITU_BS_NRT_OBSERVATIONS_013_034')
  ncatt_put(ncout,0,'licence'     ,'http://marine.copernicus.eu/services-portfolio/service-commitments-and-licence/'  ) #please do not change
  ncatt_put(ncout,0,'product'     ,'BS_OMI_OC')
  ncatt_put(ncout,0,'dataset'     ,'bs_omi_oc')
  ncatt_put(ncout,0,'quality_information_document','http://marine.copernicus.eu/documents/QUID/CMEMS-OMI-QUID-BS-OC-TSERIES.pdf')
  ncatt_put(ncout,0,'product_user_manual'         ,'http://marine.copernicus.eu/documents/PUM/CMEMS-OMI-PUM-OC.pdf')
  ncatt_put(ncout,0,'area'                        ,'BS')
  ncatt_put(ncout,0,'comment'                     ,'Period : 2010-2018')
  nc_close(ncout)
}



