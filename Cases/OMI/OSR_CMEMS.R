# Generation of the Black Sea Oxygen content files from Argo data for CMEMS OMI
# A. Capet 02/07/2018 - acapet@uliege.be

# Load package and function sources
source("../../Utils/ArgoLoad.R")
# Overwrites default ArgoSelect function
source("ArgoSelect_CMEMS.R")

datacase <- "CMEMS_MONTHLY" #"CORIOLIS" #"CMEMS_HISTORY" #  

################################
# Criterium for Argo selection #
################################

source(paste0( "OmiDataLoad_", datacase,".R" ) )

########################
# Function definitions #
########################
# To add new function, use :  
source("ArgoVertFunctions_OMI.R")
source("../../Utils/ArgoCompleteBIS.R")
source("../../Utils/ArgoExtractForCastedDF.R")

########################
# Complete (by level)  #
#      and / or        #
# Extract (by profile) #
########################

ddi  <- dcast(fuldf,qc+depth+juld+lon+lat+Platform+day+month+year~variable,fun.aggregate = mean, na.rm=T)
ddin <- subset(ddi,!is.na(DOXY))
ddim <- melt(ddin,id.vars = c("qc","depth","juld","lon","lat","Platform","day","month","year"))

flistl <- list("InSituDens")
fdf   <- ArgoCompleteBIS(ddim, flistl)

ddi<-dcast(fdf,qc+depth+juld+lon+lat+Platform+day+month+year~variable,fun.aggregate = mean, na.rm=T)
ddin<-subset(ddi,!is.na(DOXY) & !is.na(RHO))
ddim<-melt(ddin,id.vars = c("qc","depth","juld","lon","lat","Platform","day","month","year"))

flist   <- list("Z20","R20","VOC","CCC")

finaldf <- ArgoExtractForCastDF(ddin, flist )

varsdf<-rbind(
  data.frame(row.names="VOC", ymin=5,    ymax= 35,   yminl=8,    ymaxl= 35,   Fact= 1/1000 , Unit="mol m-2", YLAB= paste("Oxygen Inventory "  ,sep="")),       
  data.frame(row.names="R20", ymin=16.2, ymax= 14.6, yminl=16.3, ymaxl= 15.0, Fact= 1      , Unit="kg m-3" , YLAB= paste("Penetration density",sep="")),       
  data.frame(row.names="Z20", ymin=160,  ymax= 20,   yminl=170,  ymaxl= 70,   Fact=1       , Unit="m",       YLAB= paste("Penetration depth"  ,sep=""))      
)

finaldffrop<-dcast(finaldf,qc+juld+lon+lat+Platform+day+month+year~variable, fun.aggregate = mean, na.rm=TRUE)

######################
# OMI netcdf ouptuts #
######################

# save(finaldffrop,file = 'Lastfinaldrop.Rdata') # load(file = 'Lastfinaldrop.Rdata')
# save(finaldffrop,file = paste0('Lastfinaldrop',datacase,'.Rdata')) # load(file = paste0('Lastfinaldrop',datacase,'.Rdata'))
source("NetcdfOutput_OMI_MultiPlatform.R")
source("NetcdfOutput_OMI_MultiPlatformMonthly.R")

finaldffrop$Platform<-gsub(" ", "", finaldffrop$Platform, fixed = TRUE)
NetcdfOutput_OMI_MP(finaldffrop, 'Z20')
NetcdfOutput_OMI_MP(finaldffrop, 'R20')
NetcdfOutput_OMI_MP(finaldffrop, 'VOC')

# Procedure For OMI preparation stops here.

break()

#########
# PLOTS #
#########
finaldffrop<-subset(finaldffrop,year>=2011)

vavar<-"Z20"
ZPLOTz <-
  ggplot(finaldffrop,
         aes(x=as.Date(juld,origin = as.Date("01 jan 1950", "%d %b %Y")),y=Z20))+
  geom_point(aes(color=Platform, fill=Platform), alpha=0.8)+
  geom_smooth(aes(color=Platform, fill=Platform), method="loess",level=0.99)+
  geom_smooth(color="black",method="loess",size=2,span=0.5,level=0.99)+
  coord_cartesian(ylim = as.numeric(varsdf[vavar,c("ymin","ymax")]))+scale_y_reverse()+
  ylab(paste (varsdf[vavar,"YLAB"],"-","[",varsdf[vavar,"Unit"],"]"))+
  xlab("")+
  scale_fill_discrete(name="ARGO ID")+
  scale_colour_discrete(name="ARGO ID")+
  theme_light()+
  scale_x_date()+
  theme(legend.position = "none")
ZPLOTz

vavar<-"R20"
RPLOTz <-
  ggplot(finaldffrop,
         aes(x=as.Date(juld,origin = as.Date("01 jan 1950", "%d %b %Y")),y=R20))+
  geom_point(aes(color=Platform, fill=Platform), alpha=0.8)+
  geom_smooth(aes(color=Platform, fill=Platform), method="loess",level=0.99)+
  geom_smooth(color="black",method="loess",size=2,span=1,level=0.99)+
  scale_y_reverse(limits=c(16.25,14.75))+
  ylab(paste (varsdf[vavar,"YLAB"],"-","[",varsdf[vavar,"Unit"],"]"))+
  xlab("")+
  scale_fill_discrete(name="ARGO ID")+
  scale_colour_discrete(name="ARGO ID")+
  theme_light()+
  theme(legend.position = "none")
RPLOTz

vavar<-"VOC"
VPLOTz <-
  ggplot(finaldffrop,
         aes(x=as.Date(juld,origin = as.Date("01 jan 1950", "%d %b %Y")),y=VOC))+
  geom_point(aes(color=Platform, fill=Platform), alpha=0.8)+
  geom_smooth(aes(color=Platform, fill=Platform), method="loess",level=0.99)+
  geom_smooth(color="black",method="loess",size=2,span=1,level=0.99)+
  ylab(paste (varsdf[vavar,"YLAB"],"-","[",varsdf[vavar,"Unit"],"]"))+
  xlab("")+
  scale_fill_discrete(name="ARGO ID")+
  scale_colour_discrete(name="ARGO ID")+
  theme_light()+
  theme(legend.position = "none")
VPLOTz

require(gridExtra)
pgs<-grid.arrange(ZPLOTz,RPLOTz,VPLOTz, ncol = 1, nrow = 3)

source("~/Desktop/PAPERS_UNDER_WORK/OSR3/scripts/g_legend.R")
plegend <- 
  ggplot(finaldffrop,
         aes(x=as.Date(juld,origin = as.Date("01 jan 1950", "%d %b %Y")),y=VOC))+
  geom_point(aes(color=Platform, fill=Platform), alpha=0.8)+
  geom_smooth(aes(color=Platform, fill=Platform), method="loess",level=0.99)+
  geom_smooth(color="black",method="loess",size=2,span=1,level=0.99)+
  scale_y_reverse()+
  ylab(paste ("Depth -","[",varsdf[vavar,"Unit"],"]"))+
  xlab("")+
  scale_fill_discrete(name="ARGO ID")+
  scale_colour_discrete(name="ARGO ID")+
  theme_light()

plegend<-g_legend(plegend)

pdf(file="~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/CMEMS_F2b.pdf", 
    width  = 10,
    height = 15,#8.27 / (1.68-0.2),
    family = "Times")
grid.arrange(pgs,plegend, ncol = 2, nrow = 1,widths=c(4,1))
dev.off()

pdf(file="~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/CMEMS_F2bOldStyle.pdf", 
    width  = 12,
    height = 12 / (1.68-0.2),
    family = "Times")
grid.arrange(pgs,plegend, ncol = 2, nrow = 1,widths=c(4,1))
dev.off()


###########
# FIG2    #
###########
library(ggmap)
blackseamap <- get_map(location = c(27.5,41,42,47),source="google",color = "bw")

  p1<-
    ggmap(blackseamap)+
    geom_path(data=finaldffrop, aes(x=lon, y=lat,color=Platform))+
    geom_point(data=finaldffrop, aes(x=lon, y=lat,color=Platform))+
    facet_wrap(~year)+theme_light()+
    scale_colour_discrete(name="ARGO ID")
  p2<-
    ggplot(finaldffrop, aes(x=as.Date(juld,origin = as.Date("01 jan 1950", "%d %b %Y")),
                            y=Platform,color=Platform))+
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")#+
  #  scale_y_discrete(name="ARGO ID",labels=ARGOLAB)
  
#  pdf(paste("~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/Map.pdf",sep=""),width=15, height=8)
  pdf(paste("~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/CMEMS_Map.pdf",sep=""),width=8, height=15)
  grid.arrange(p2,p1,ncol=1)#, heights=c(6,10))
  dev.off()

########
# FIG3 #
########

##########################
## Getting CTD diagnostics
##########################
  
  DDir<-"/home/arthur/diva/BlackSea4diva/"
  fi<-paste(DDir,"data/","VOC","/detrending/","VOC",".dat",sep="")
  CTDdata<-read.table(fi,col.names=c("lat","lon","Value","Month","Year","Weight"))
  
  for (VAR in c('CCO','VOC')){
    fi<-paste(DDir,"data/",VAR,"/detrending/",VAR,".dat",sep="")
    totdata<-read.table(fi,col.names=c("lat","lon","Value","Month","Year","Weight"))
    CTDdata[,VAR]<-totdata[,"Value"]
    CTDdata[,"Set"]<-"CTD"
  }

  yearmin<-1900
  
  totdata<-subset(totdata,Year>yearmin)   #### JUSTIFY THIS ####
  #globalaverage<-mean(totdata$Value)
  
  dfa<-finaldffrop[,c("CCC","VOC","year")]
  colnames(dfa)<-c("CCC","VOC","Year")
  dfc<-CTDdata[,c("CCO","VOC","Year")]
  dfc$VOC<-dfc$VOC/1000
  colnames(dfc)<-c("CCC","VOC","Year")
  
  df<-rbind(dfa,dfc)
  colnames(df)<-c("CCO","VOC","Year")
  
  CCCconvert<- -1014*4079/1e6 # from mK to 10^6 J/m²
  df$CCO<-df$CCO*CCCconvert
  
  df<-subset(df,!is.na(VOC))
  ################################################
  
  YBREAKS<-c(1955,1975,1985,1998,2016)
  
  for (YB in 1:(length(YBREAKS)-1)){
    CTDdata[CTDdata[,"Year"]>YBREAKS[YB] & CTDdata[,"Year"]<=YBREAKS[YB+1],"Period"] <- paste(YBREAKS[YB],"-",YBREAKS[YB+1])
    df[df[,"Year"]>YBREAKS[YB] & df[,"Year"]<=YBREAKS[YB+1],"Period"] <- paste(YBREAKS[YB],"-",YBREAKS[YB+1])
  }
  
  ######
  
  MINYEAR <- YBREAKS[1]
  OXLIM    <- c(4,42)
  VOCLIM    <- c(0,40)
  CCCLIM    <- c(1,100)
  OXLIM    <- c(7,42)
  CCCLIM    <- c(0,100)*CCCconvert
  XLAB    <- "CIL Cold Content - [10^6 J/m^2]"
  YLAB    <- "Oxygen Inventory - [mol/m^2]"
  PAL<- "Set2"# "Spectral" #"Pastel1"#, "Spectral"
  
  CTDloc<-subset(df, Year>MINYEAR & Year<=2016)
  CTDloc_2016<-subset(df, Year>=2016 & Year<2017. )
  CTDloc_2016$Period<-"2016"
  CTDloc_2017<-subset(df, Year>=2017 )
  CTDloc_2017$Period<-"2017"
  
  #################################
  
  METH = 'loess'
  SPAN = 0.75 # SPAN is used for loess, ignored otherwise
  central   <-
    ggplot(CTDloc,aes(x=CCO,y=VOC,color=Period, fill = Period))+
    geom_point(alpha=0.4)+
    geom_smooth(data=CTDloc_2017,size=2, level=0.99, linetype=1, alpha=0.7,method=METH,span=SPAN)+
  #  geom_smooth(data=CTDloc_2016,size=2, level=0.99, linetype=1, alpha=0.7,method=METH,span=SPAN)+
    geom_smooth(size=2, level=0.99, linetype=1, alpha=0.4,method=METH,span=SPAN)+
   # geom_point(data=CTDloc_2016,pch=21,color='black')+
    geom_point(data=CTDloc_2017,pch=21,color='black')+
   # geom_point(data=CTDloc_2016,pch=21,color='black')+
    ylim(OXLIM)+xlim(CCCLIM)+
    scale_color_brewer(palette = PAL, name="Period")+
    scale_fill_brewer(palette = PAL, name="Period")+
    theme_light()+
    theme(legend.text =element_text(size=13),
          legend.title =element_text(size=13))
  
  central
  
  leg<-g_legend(central)
  central   <- central  +theme(legend.position = "none",panel.grid = element_line(colour = "lightgrey"))
    
  lowband   <- ggplot(CTDloc, aes(x=CCO, fill = Period)) +
    geom_density(alpha = 0.5)+
    geom_density(data=CTDloc_2017,alpha = 0.5)+
    #geom_density(data=CTDloc_2016,alpha = 0.5)+
    #  geom_density(data=CTDloc_2016,color='black',alpha = 0.5)+
    scale_fill_brewer(palette = PAL)+ xlim(CCCLIM)+
    xlab(XLAB)+ theme_linedraw() +
    theme_light()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  leftband  <- ggplot(CTDloc, aes(x=VOC, fill = Period)) +
    geom_density(alpha = 0.5)+
    geom_density(data=CTDloc_2017,alpha = 0.5)+
    geom_density(data=CTDloc_2016,alpha = 0.5)+
    scale_fill_brewer(palette = PAL)+xlab(YLAB)+
    theme_light()+
    xlim(OXLIM)+ 
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    coord_flip()
  
  empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), 
                                                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                                      panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), 
                                                                      axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), 
                                                                      axis.ticks = element_blank())
  
  grid.arrange(leftband, central,  leg, lowband,
               ncol = 2, nrow = 2,
               widths = c(1,4), heights = c(4, 1))
  
  
  pdf(file="~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/CMEMS_OSR3_F3b.pdf", 
      width  = 8.27,
      height = 8.27 /1.68,
      family = "Times")
  grid.arrange(leftband, central,  leg, lowband,
               ncol = 2, nrow = 2,
               widths = c(1,5), heights = c(5, 2))
  dev.off()

########
# FIG4 #
########
  ll<-lapply(rownames(varsdf), function(VAR){
  DDir<-"/home/arthur/diva/BlackSea4diva/"
  fi<-paste(DDir,"data/",VAR,"/detrending/",VAR,".dat",sep="")
  totdata<-read.table(fi,col.names=c("lat","lon","Value","Month","Year","Weight"))
  totdata[,"Value"]<-totdata[,"Value"]*varsdf[VAR,"Fact"]
  totdata<-subset(totdata,Year>yearmin)   #### JUSTIFY THIS ####
  globalaverage<-mean(totdata$Value)
  
  # LOAD INTERANNUAL TREND
  fi<-paste(DDir,"results/",VAR,"/detrend_L0pt8/trends2_00.dat",sep="")
  #  fi2<-paste(DDir,"results/",VAR,"/detrend_L0pt8/trends2_err.dat",sep="")
  
  trendi<-read.table(fi,colClasses=c("NULL",NA),col.names=c("bid","Trend"))*varsdf[VAR,"Fact"]
  # trendi<-read.table(fi2,colClasses=c("NULL",NA,NA),col.names=c("bid","Trend","err"))*varsdf[VAR,"Fact"]
  
  trendi<-cbind(Year=seq(1923,2013),trendi)
  #   trendi2<-cbind(Year=seq(1923,2005),trendi2)
  # get num of data 
  ni<-data.frame(table(totdata$Year))
  colnames(ni)<-c("Year","Nprofs")
  trendi<-merge(trendi,ni,all.x=TRUE) #
  trendi[,"Trend"]<-trendi[,"Trend"]+globalaverage
  trendi[,"VAR"]<-VAR
  return(trendi)
  })
  
  trendf<-do.call(rbind,ll)
  
  trendf[is.na(trendf$Nprofs),"Nprofs"]<-0
  
  lm_report = function(df){
    m = lm(Trend ~ Year, df);
    eq <- paste(format(coef(m)[2]*10, digits = 2),
                Unit,"/decades  R2a=", format(summary(m)$adj.r.squared, digits = 3),
                " p=", format(summary(m)$coefficients["Year",4], digits = 3)
    )
    as.character(as.expression(eq));                 
  }
  
  
  finaldffroploc      <- finaldffrop
  finaldffroploc$Set  <- "Argo"
  finaldffroploc$time <- as.Date(finaldffroploc$juld,origin = as.Date("01 jan 1950", "%d %b %Y"))
  finaldffroploc$time <-2010+julian(finaldffroploc$time, origin = as.Date("2010-01-01"))/365.25
  
  Nprofmin <- 5
  yearmin <- 1955
  yearmax <- 2006
  
  plist<-lapply(c("Z20","R20","VOC"), function(vavar) {
    ggplot(subset(trendf,Nprofs>Nprofmin&Year>yearmin&VAR==vavar&Year<yearmax)
           ,aes(x=Year,y=Trend))+
      geom_point(size=3,aes(shape="CTD"))+
      geom_line(linetype=2)+
      stat_smooth(linetype=0)+
      theme_bw()+#theme(legend.position="bottom")+
      geom_smooth(method = "lm", se=FALSE, color="red", formula = y ~ x)+
      geom_smooth(data = subset(finaldffroploc,time>2011), aes_string(x="time",y=vavar, fill="Set", color="Set"), method='loess', span=0.75, level=0.99)+
      ylim(as.numeric(varsdf[vavar,c("ymin","ymax")]))+
      ylab(paste (varsdf[vavar,"YLAB"],"-","[",varsdf[vavar,"Unit"],"]"))+
      scale_x_continuous(breaks=seq(from=1960,to=2020,by=10))
  })
  
  legend<-g_legend(plist[[1]])
  
  plist[[1]]<-plist[[1]]+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  plist[[2]]<-plist[[2]]+ theme(axis.text.x = element_blank(),axis.title.x = element_blank())
  plist<-lapply(plist,function(p) p+ theme(legend.position="none")) # ,plot.margin=unit(c(.1,.1,.1,.1), "cm")
  
  basedir<-"~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/"
  pdf(paste(basedir,VAR,yearmin, "CMEMS_F4.pdf",sep=""),width=10,heigh=6*1.68)
  do.call(grid.arrange,plist)
  dev.off()
  