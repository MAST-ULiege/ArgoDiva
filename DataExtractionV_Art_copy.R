# A. Capet : Argo data extraction : Examples
  
# SECTION ONE : EXTRACTING BASIC DATA 

library(ncdf4)
#More complex libraries for plotting data
library(ggplot2)
library(plyr)
#for melt function
library(reshape2)
  
filename <- "6901866_Mprof.nc"

#Opening the file in a open-only mode
ncfile <- nc_open(filename, write = FALSE, verbose = TRUE, suppress_dimvals = FALSE)

#Dimensions
N_PROF <- ncol(ncvar_get(ncfile,"PRES"))
N_LEVELS <- nrow(ncvar_get(ncfile,"PRES"))

juld <- ncvar_get(ncfile,"JULD")
pres <- ncvar_get(ncfile,"PRES")

##########################
# one function to extract variable and return a dataframe

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
  return(d=d)
}  
  
TEMPL<-ExtractVar("TEMP")
DOXYL<-ExtractVar("DOXY")
CDOML<-ExtractVar("CDOM")
PSALL<-ExtractVar("PSAL")

# SECTION THREE : PLOTS OF DATA ACCORDING QC VALUES -------------------------------

# Faire des points pour les valeur de température, un panneau par valeur de QC (facet_wrap)
ggplot(TEMPL,aes(x=value,y=-depth,color=aprofile))+
  geom_point()+facet_wrap(~qc)+theme_bw()+ylim(c(-200,0))

# diviser les profiles en 10 groupes egaux sur base du jour, tracer un "smoother" plutot que les points
ggplot(DOXYL,aes(y=value,x=-depth,color=cut(juld,breaks = 10)))+
  geom_smooth()+facet_wrap(~qc)+theme_bw()+coord_flip()

# diviser les profils en 10 groupes sur base du numéro de profil pour différentes variables en un coup
ggplot(subset(rbind(DOXYL,TEMPL),qc!=0),aes(y=value,x=-alevel,color=cut(aprofile,breaks = 10)))+
  geom_point()+facet_grid(qc~variable,scales = "free")+theme_bw()+coord_flip()

# ne garder que les QC = 3
ggplot(subset(CDOML,qc==3) ,aes(x=value,y=-depth,color=juld))+
  geom_point()+theme_bw()

################################
# Combining several variables: 
DOXYLc<-subset(DOXYL,qc=1,select=c("depth","juld","value"))
colnames(DOXYLc)[which(colnames(DOXYLc)=="value")]<-"DOXY"

TEMPc<-subset(TEMPL,qc=1,select=c("depth","juld","value"))
colnames(TEMPc)[which(colnames(TEMPc)=="value")]<-"TEMP"

PSALLc<-subset(PSALL,qc=1,select=c("depth","juld","value"))
colnames(PSALLc)[which(colnames(PSALLc)=="value")]<-"PSAL"

bi<-merge(DOXYLc,TEMPc,by=c("depth","juld"))
bi<-merge(bi,PSALLc,by=c("depth","juld"))

# T, S diagrams, with color for oxygen, split in 9 panels with different periods
# (with equally distributed number of profiles) 
ggplot(bi,aes(y=PSAL, x=TEMP, color=DOXY))+
  geom_point()+facet_wrap(~cut(juld,9))+theme_bw()

################################
# Executing a function on each profile

# This command separates "bi"into subdataframes, based on the values of "juld"
# computes the new variables defined after "summarize" by applying computations
# on the each resulting subdataframe.
# Then it merges those back into a large dataframe
maximums <- ddply(bi,~juld,summarize, 
                  maxTem=max(TEMP),
                  maxDoxy=max(DOXY), 
                  depthofmaxDoxy = depth[which.max(DOXY)])

ggplot(melt(maximums,id.vars = "juld"),aes(x=juld,y=value, color=variable))+
  geom_path()+
  facet_wrap(~variable,ncol=1,scales = "free_y")

ggplot(maximums,aes(x=maxTem,y=maxDoxy, color=juld))+
  geom_point()
