NetcdfOutput_OMI <- function(df,varname){
  
  # -- Path here to save the netcdf file
  path <- paste0('./Output_',Sys.Date())
  
  # -- OMI file name : let's build the standard name of your future file
  area            <- 'bs'             # or med or nws, ...
  short_omi_name  <- 'oc'             # 
  first_occurence <- 'FIRSTOCCURENCE' # yyyymmdd
  prod_time       <- paste0('P',gsub("-",'',as.character(Sys.Date())))      # P + yourdate yyyymmdd
  other_info      <- ''               # ohc case ...
  
  ncfname <- paste0(area,'_omi_',short_omi_name,'_area_averaged_',first_occurence,other_info,'.nc')
  
  # -- TIME info
  TIMEdf <- data.frame(origin   = as.Date('1950-01-01 00:00:00'),
                       calendar = 'gregorian',
                       units    = 'days since 1950-01-01 00:00:00',
                       name     = 'time')

  # -- VAR infos
  VOCdf   <- data.frame(long_name     = "Oxygen Content in a water column",
                       units         = 'mol/m2',
                       name          = 'oxygen_content',
                       standard_name = 'oxygen_content',
                       FillValue     = -9999.0,
                       Comment       = 'Vertically integrated oxygen concetration')
  rownames(VOCdf)<-"VOC"
  
  Z20df   <- data.frame(long_name     = "Oxygen Penetration Depth",
                       units         = 'm',
                       name          = 'oxycline_depth',
                       standard_name = 'oxycline_depth',
                       FillValue     = -9999.0,
                       Comment       = 'Depth at which [O2]<20µM')
  rownames(Z20df)<-"Z20"
  
  R20df   <- data.frame(long_name     = "Oxygen Penetration Density",
                        units         = 'kg/m3',
                        name          = 'oxycline_density',
                        standard_name = 'oxycline_density',
                        FillValue     = -9999.0,
                        Comment       = 'Potential density anomaly at which [O2]<20µM')
  rownames(R20df)<-"R20"
  
  vardf<-rbind(OCdf,ODNdf,ODPdf)
  
  # -- Create the dimensions of the netcdf file
  timedim <- ncdim_def(TIMEdf$name,TIMEdf$units, unlim = TRUE)

  # -- Create 
  for(v in rownames(vardf)){
    ncvar_def(name     = vardf[v,'name'],
              units    = vardf[v,'units'],
              dim      = list(timedim),
              missval  = vardf[v,'FillValue'], 
              longname = vardf[v,'long_name'])    
  }

  # -- Create netcdf file  
  ncout   <- nc_create(ncfname,list(AR_def,MO_def,AT_def,DI_def,CO_def),force_v4=T)
  
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
    ncatt_put(ncout,vardf[v,'name'], "standard_name",vardf[v,'standard_name'])
    ncatt_put(ncout,vardf[v,'name'], "comment",vardf[v,'comment'])
  }
  

#  Global Attributes
  ncatt_put(ncout,0,'title')           <- 'Black Sea Oxygen Inventory' # to be consistent with PIT
  ncatt_put(ncout,0,'institution')                  <- 'ULiege'
  ncatt_put(ncout,0,'references')                   <- 'http://marine.copernicus.eu'   #please do not change
  ncatt_put(ncout,0,'Conventions')                  <- 'CF-1.7'
  ncatt_put(ncout,0,'credit')                       <- 'E.U. Copernicus Marine Service Information (CMEMS)' #please do not change
  ncatt_put(ncout,0,'contact')                      <- 'servicedesk.cmems@mercator-ocean.eu'                #please do not change
  ncatt_put(ncout,0,'source')                       <- 'INSITU_BS_NRT_OBSERVATIONS_013_034'
  ncatt_put(ncout,0,'licence')                      <- 'http://marine.copernicus.eu/services-portfolio/service-commitments-and-licence/'   #please do not change
  ncatt_put(ncout,0,'product')                      <- 'GLOBAL_OMI_OHC_area_averaged_anomalies'
  ncatt_put(ncout,0,'dataset')                      <- 'global_omi_ohc_area_averaged_anomalies'
  ncatt_put(ncout,0,'quality_information_document') <- 'http://marine.copernicus.eu/documents/QUID/CMEMS-OMI-QUID-GLO-OHC-TSERIES.pdf'
  ncatt_put(ncout,0,'product_user_manual')		      <- 'http://marine.copernicus.eu/documents/PUM/CMEMS-OMI-PUM-OHC.pdf'
  ncatt_put(ncout,0,'area')	 	         		          <- 'BS'
  ncatt_put(ncout,0,'comment')			                <- 'Period : 2010-2018'
}