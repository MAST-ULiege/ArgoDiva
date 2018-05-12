# Case Test example for the use of the ArgoDiva toolbox
# A. Capet 02/02/2018

# Load package and function sources
source("ArgoLoad.R")

################################
# Criterium for Argo selection #
################################

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2010-2014/"
selectCriterium0<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
fulldf0 <- ArgoSelect(selectCriterium0)

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2014-2016/"
selectCriterium1<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
fulldf1 <- ArgoSelect(selectCriterium1)

dataDir   = "/home/arthur/Desktop/PAPERS_UNDER_WORK/OSR3/datas/2017/"
selectCriterium2<-list(
  dataDir   = dataDir,
  ## List of individual Argo netcdf files. 
  # filenames = c("6901866_Mprof.nc"),
  # filenames = c("argo-profiles-5902291.nc","argo-profiles-6901960.nc"),
  ########################################
  ## To get all files in the directory use 
  filenames = list.files(dataDir,"argo-profiles-.*\\.nc"),#[1:5],
  #filenames = Sys.glob(paste0(dataDir,"argo-profiles-*.nc")),
  ## List of variable to be extracted from Argo  ( should correspond to Argo variable names)
  varList      = c("TEMP","PSAL","DOXY"), 
  # This flag states wether we retain only records for which all variables are presents
  presentInAll = TRUE
)
fulldf2 <- ArgoSelect(selectCriterium2)

fuldf<-rbind(fulldf0,fulldf1,fulldf2)
########################
# Function definitions #
########################

# To know about defined function you could use , look at :
source("ArgoVertFunctions.R")

# To add new function, use :  (using same structure as in the previous file)
source("ArgoVertFunctions_USER.R")

########################
# Complete (by level)  #
#      and / or        #
# Extract (by profile) #
########################

ddi  <- dcast(fuldf,qc+depth+juld+lon+lat+Platform+day+month+year~variable,fun.aggregate = mean, na.rm=T)
ddin <- subset(ddi,!is.na(DOXY))
ddim <- melt(ddin,id.vars = c("qc","depth","juld","lon","lat","Platform","day","month","year"))

flistl <- list("InSituDens")
source("ArgoCompleteBIS.R")
fdf   <- ArgoCompleteBIS(ddim, flistl)

#ddi<-subset(fdf,  Platform=="6901866 ")
#require(ggplot2)
#ggplot(subset(ddi,aprofile>50&aprofile<60),aes(x=value,y=depth))+
#  geom_point()+scale_y_reverse()+facet_grid(.~variable, scales = "free")+ylim(c(32,30))

#ddi2<-dcast(ddi,qc+alevel+depth+aprofile+juld+lon+lat+Platform+day+month+year~variable)
ddi<-dcast(fdf,qc+depth+juld+lon+lat+Platform+day+month+year~variable,fun.aggregate = mean, na.rm=T)
ddin<-subset(ddi,!is.na(DOXY) & !is.na(RHO))
ddim<-melt(ddin,id.vars = c("qc","depth","juld","lon","lat","Platform","day","month","year"))

flist <- list("Z20","R20","VOC","CCC")
source("ArgoExtractForCastedDF.R")
finaldf   <- ArgoExtractForCastDF(ddin, flist )

varsdf<-rbind(
  data.frame(row.names="VOC", ymin=5, ymax= 35,  yminl=8, ymaxl= 35, Fact= 1/1000 , Unit="mol/m²", YLAB= paste("Oxygen Inventory ",sep="")),       
  data.frame(row.names="R20", ymin=16.2, ymax= 14.6, yminl=16.3, ymaxl= 15.0 , Fact= 1      , Unit="kg/m³" , YLAB= paste("sigma 20",sep="")),       
  data.frame(row.names="Z20", ymin=160,ymax= 20, yminl=170,ymaxl= 70, Fact=1, Unit="m", YLAB= paste("z 20",sep=""))      
)

finaldffrop<-dcast(finaldf,qc+juld+lon+lat+Platform+day+month+year~variable, fun.aggregate = mean, na.rm=TRUE)


##########
# PLOTS
##########
finaldffrop<-subset(finaldffrop,year>=2011)

vavar<-"Z20"
ZPLOTz <-
  ggplot(finaldffrop,
         aes(x=as.Date(juld,origin = as.Date("01 jan 1950", "%d %b %Y")),y=Z20))+
  geom_point(aes(color=Platform, fill=Platform), alpha=0.8)+
  geom_smooth(aes(color=Platform, fill=Platform), method="loess",level=0.99)+
  geom_smooth(color="black",method="loess",size=2,span=0.5,level=0.99)+
  coord_cartesian(ylim = as.numeric(varsdf[vavar,c("ymin","ymax")]))+scale_y_reverse()+
  ylab(paste ("Depth -","[",varsdf[vavar,"Unit"],"]"))+
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
  ylab(paste ("Depth -","[",varsdf[vavar,"Unit"],"]"))+
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
  #scale_y_reverse()+
  ylab(paste ("Depth -","[",varsdf[vavar,"Unit"],"]"))+
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


pdf(file="~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/F2b.pdf", 
    width  = 10,
    height = 15,#8.27 / (1.68-0.2),
    family = "Times")
grid.arrange(pgs,plegend, ncol = 2, nrow = 1,widths=c(4,1))
dev.off()

pdf(file="~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/F2bOldStyle.pdf", 
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
  pdf(paste("~/Desktop/PAPERS_UNDER_WORK/OSR3/figure/Map.pdf",sep=""),width=8, height=15)
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
  globalaverage<-mean(totdata$Value)
  
  # Here You Are #
  # -> Get CCC from current Argo 
  # -> Adapt time issues
  
  dfa<-ARGO[,c("CCC","VOC","time")]
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
  
  CTDloc<-subset(df, Year>MINYEAR & Year<2016)
  CTDloc_2016<-subset(df, Year>2016 & Year<=2017. )
  CTDloc_2016$Period<-"2016"
  CTDloc_2017<-subset(df, Year>2017 )
  CTDloc_2017$Period<-"2017"
  
  #################################
  
  METH = 'loess'
  SPAN = 0.75 # SPAN is used for loess, ignored otherwise
  central   <-
    ggplot(CTDloc,aes(x=CCO,y=VOC,color=Period, fill = Period))+
    geom_point(alpha=0.5)+
    geom_smooth(size=2, level=0.99, linetype=1, alpha=0.6,method=METH,span=SPAN)+
    ylim(OXLIM)+xlim(CCCLIM)+
    scale_color_brewer(palette = PAL, name="Period")+
    scale_fill_brewer(palette = PAL, name="Period")+
    theme_light()+
    theme(legend.text =element_text(size=13),
          legend.title =element_text(size=13))
  leg<-g_legend(central)
  
  central   <-
    ggplot(CTDloc,aes(x=CCO,y=VOC,color=Period, fill = Period))+
    geom_point(alpha=0.5)+
    geom_point(data=CTDloc_2016,color='black')+
    geom_point(data=CTDloc_2017,color='red')+
    geom_smooth(size=2, level=0.99, alpha=0.6,method=METH,span=SPAN)+
    ylim(OXLIM)+xlim(CCCLIM)+
    scale_color_brewer(palette = PAL, name="")+scale_fill_brewer(palette = PAL, name="")+
    xlab(XLAB)+ylab(YLAB)+
    theme_light()+
    theme(legend.position = "none",panel.grid = element_line(colour = "lightgrey"))
  
  
  lowband   <- ggplot(CTDloc, aes(x=CCO, fill = Period)) +
    geom_density(alpha = 0.5)+
    #  geom_density(data=CTDloc_2016,color='black',alpha = 0.5)+
    scale_fill_brewer(palette = PAL)+ xlim(CCCLIM)+
    xlab(XLAB)+ theme_linedraw() +
    theme_light()+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  
  leftband  <- ggplot(CTDloc, aes(x=VOC, fill = Period)) + geom_density(alpha = 0.5)+
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
  
  pdf(paste("../figure/OSR3_VOC_CCC_allin_NN_linear_2107.pdf",sep=""), height = 7,width=7*1.68)
  grid.arrange(leftband, central,  leg, lowband,
               ncol = 2, nrow = 2,
               widths = c(1,5), heights = c(5, 2))
  dev.off()
  
  
  pdf(file="../figure/OSR3_Sec3.5_F3.pdf", 
      width  = 8.27,
      height = 8.27 /1.68,
      family = "Times")
  grid.arrange(leftband, central,  leg, lowband,
               ncol = 2, nrow = 2,
               widths = c(1,5), heights = c(5, 2))
  dev.off()
  
  
  postscript(file="../figure/OSR3_Sec3.5_F3.eps", 
             width  = 8.27,
             height = 8.27 /1.68,
             family = "Times")
  grid.arrange(leftband, central,  leg, lowband,
               ncol = 2, nrow = 2,
               widths = c(1,5), heights = c(5, 2))
  dev.off()
  
  
  
  


#Visualization of distribution density for different period (showing different approach)
ArgoDisplay(fdf,selectCriterium,flist)

##################################
# Parameters for diva input file #
##################################



