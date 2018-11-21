###########################################
# Functions for "levels", ie. depths

InSituDens <- function(llevel){
  require(gsw)
  tem <- subset(llevel,variable=="TEMP", value)
  psal <- subset(llevel,variable=="PSAL", value)

  psal <- gsw_SA_from_SP(psal,unique(llevel$depth),unique(llevel$lon),unique(llevel$lat))
  rho_anomaly <- gsw_sigma0(psal,tem)
  
  out <- data.frame( value    = rho_anomaly, qc = min(llevel$qc),  variable = 'RHO')
  return(out)
}

###########################################
# Functions for profiles


# Needed for PAR Works

grad_att <- function(PRESr,k,PAR0){
  PAR<-PRESr*0
  kint<-(c(k[1],k)+c(k,k[length(k)]))/2
  PAR[1]<-PAR0*exp(-kint[1]*PRESr[1])
  for (i in 2:length(PAR)){
    PAR[i]<-PAR[i-1]*exp(-((kint[i]+kint[i])/2    )*(PRESr[i]-PRESr[i-1]))
  }
  return(PAR)
}

Att_2band <- function(PRES,PAR0,ks0,kl0,depthP0=0) {
  ##  "CHL_2BANDS" ##
  # CHLA specific attenuation + two bandwidths
  if (depthP0 !=0){
    PAR0  = PAR0/(parts*exp(-ks0*depthP0)+(1-parts)*exp(-kl0*depthP0)  ) 
  }
  PAR0s <- PAR0*parts
  PAR0l <- PAR0*(1-parts)
  ks    <- rep(ks0,length(PRES))
  kl    <- rep(kl0,length(PRES))
  PARs  <- grad_att(PRES,ks,PAR0s)
  PARl  <- grad_att(PRES,kl,PAR0l)
  PAR<-PARs+PARl
  return(PAR)
}





p2b_start = c(PAR0 = 1500 ,  ks0 = 0.04 , kl0=0.001) 
p2b_lower = c(PAR0 = 0    ,  ks0 = 0    , kl0=0    )
p2b_upper = c(PAR0 = 2500 ,  ks0 = 3   , kl0=10)

Out_Default <- data.frame("PAR0"=NA,
                          "ks0"=NA,
                          "kl0"=NA)

parts = 0.37

## Too Slow
# require(FME)
# f2bandFME <- function(p, profiledata){
#   with(as.list(p),{
#     MOD <- Att_2band(profiledata$depth,PAR0,ks0,kl0)
#     OBS <- profiledata$PAR
# #    return(modCost(model=MOD, obs = OBS))
#     return(MOD-OBS)
#     
#   })
# }
# 
# PAR2bands_FME <- function(profile){
#   # Intended for Casted DF
#   Out_Default <- data.frame("PAR0"=NA,
#                             "parts"=NA,
#                             "kl0"=NA)
#   if (FilterForPAR(profile)){
#     out <-try_default({
#       
#       m   <- modFit(f = f2bandFME,
#                     p = p2b_start,
#                     lower = p2b_lower , 
#                     upper = p2b_upper , 
#                     method = "Pseudo" ,
#                     profiledata = profile )
#         
#       return(as.data.frame(t(coef(m))))},
#       Out_Default)
#   }else{
#     out <- Out_Default
#   }
#   return(out)
# }
#####################################################
# NLS #

PAR2bands <- function(profile){
  # Intended for Casted DF
  if (FilterForPAR(profile)){
    out <-try_default({
      m   <-nls(PAR ~ Att_2band(depth,PAR0,ks0,kl0),
          data = profile,
          start = p2b_start, 
          lower = p2b_lower,
          upper = p2b_upper,
          algorithm = "port", 
          control = list(maxiter= 1000)
          )
      return(as.data.frame(t(coef(m))))},
    Out_Default)
  }else{
    out <- Out_Default
  }
  return(out)
}


DefOut_isoPAR <- data.frame("z10"=NA,
                            "z100"=NA)

isoPAR <- function(profile){
  # Intended for Casted DF
  if (FilterForPAR(profile)){
    out    <- try_default({


    d <- profile$depth[which( !is.na(profile$depth) & profile$PAR>0.1  )]
    P <- profile$PAR[which( !is.na(profile$depth) & profile$PAR>0.1  )]

    p1 <-interp1(d,P, 1)      
#      pcP <- pchip(profile$depth,profilePAR, seq(1,200))
      
      z10  <- interp1(P,d, p1/10) 
      z100 <- interp1(P,d, p1/100) 
      return(data.frame(z10=z10,z100=z100))},
      DefOut_isoPAR)
  }else{
    out <- DefOut_isoPAR
  }
  return(out)
}





#####################################################
R20 <- function(profile){
  # Intended for Casted DF
  if (FilterForVOX(profile)){
    out <- data.frame( variable = "R20",value=min(profile$RHO[which(profile$DOXY<20)],na.rm = T))
    }else{
    out <- data.frame( variable = "Z20", value=NA)
    }
  return(out)
}

Z20 <- function(profile){
  # Intended for Casted DF
  if (FilterForVOX(profile)){
    
  out <- data.frame(
    variable = "Z20",
    value=min(profile$depth[which(profile$DOXY<20)],na.rm = T)
  )
  if (out$value<20 |out$value>300){
    plot(profile$DOXY,profile$depth,ylim=c(500,0) )
    lines(pchip(profile$depth,profile$DOXY, seq(1,300)),seq(1,300))
    lines(c(20,20),c(400,0), col='red')
    print(out$value)
    print(profile)
    print(out$value)
  }
  }else{
    out <- data.frame(
      variable = "Z20",
      value=NA
    )
  }
  return(out)
}

# TODO limit to upper 20µM
require(pracma)
VOC <- function(profile){
  # Intended for Casted DF
  if (FilterForVOX(profile)){
    Z20<-Z20(profile)$value
    d  <- profile$depth[which(profile$depth<=Z20)]
    o  <- profile$DOXY[which(profile$depth<=Z20)]
    zp <- seq(1,Z20)
    po <- pchip(d,o, seq(1,Z20))
    po[which(zp<min(d))]<-o[which(d==min(d))]
    try_default(
  out <- data.frame(
    variable = "VOC",
    value=sum(po)/1000
    ),
  out <- data.frame(
    variable = "VOC",
    value=NA  ), 
  TRUE
  )
     if (is.na(out$value)|out$value>5e4|out$value<100){
       print(profile)
print(       min(profile$depth))
     plot(profile$DOXY,profile$depth,ylim=c(200,0) )
      lines(c(20,20),c(400,0), col='red')
      lines(pchip(d,o, seq(1,Z20)),seq(1,Z20))
      lines(po,zp,col='green')
      points(pchip(d,o, seq(1,Z20)),seq(1,Z20), pch=3)
      print(out$value)
      }
  }else{
    out <- data.frame(
      variable = "VOC",
      value=NA
    )
  }
  return(out)
}

CCC <- function(profile){
  # Intended for Casted DF
  if (FilterForVOX(profile)){
    #Z20<-Z20(profile)$value
    d  <- profile$depth#[which(profile$depth<=Z20)]
    t  <- profile$TEMP#[which(profile$depth<=Z20)]
    r  <- profile$RHO
    zp <- seq(1,200)
    pt <- pchip(d,t, zp)
    pr <- pchip(d,r, zp)
    #po[which(zp<min(d))]<-o[which(d==min(d))]
    pt <- 8.35-pt
    pt[which(pt<=0)]<-0
    pt[which(pr<=14)]<-0
    try_default(
      out <- data.frame(
        variable = "CCC",
        value=sum(pt)
      ),
      out <- data.frame(
        variable = "CCC",
        value=NA  ), 
      TRUE
    )
    if (is.na(out$value)|out$value>10000|out$value<0){
      print(profile)
      print(min(profile$depth))
      plot(profile$TEMP,profile$depth,ylim=c(200,0) )
      lines(c(8.35,8.35),c(400,0), col='red')
      lines(pchip(d,t, zp),zp,col='green')
      print(out)
    }
  }else{
    out <- data.frame(
      variable = "CCC",
      value=NA
    )
  }
  return(out)
}

###########################################
# Criterion for function evaluation
#
# Given below are conditions for exclusion 

FilterForVOX <- function(profile){
  # Gather here requirements for Oxygen diagnostic extraction
  flag=TRUE
  if (sum(!is.na(profile$DOXY))<10 | min(profile$DOXY)>20 | min(profile$depth)>15){
    flag=FALSE
  }
  return (flag)
}

FilterForCCC <- function(profile){
  # Gather here requirements for CCC diagnostic extraction
  flag=TRUE
  if (sum(!is.na(profile$TEMP)) < 10 |   # min 10 depths
          #min(profile$TEMP)>8.35   |   # CIL is defined
          #min(profile$depth)>15    |   # CIL is defined
          max(profile$depth)    < 80
      ){flag=FALSE}
  return (flag)
}

FilterForPAR <- function(profile){
  # Gather here requirements for PAR diagnostic extraction
  flag=TRUE
  if (sum(!is.na(profile$PAR)) < 10 |   # min 10 depths
      # min(profile$TEMP)>8.35   |   # CIL is defined
       min(profile$depth) > 1  #  |   # CIL is defined
      # max(profile$depth)    < 80
  ){flag=FALSE}
  return (flag)
}


