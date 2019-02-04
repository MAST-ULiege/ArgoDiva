###########################################
# Functions for "levels", ie. depths

PotDensAnom <- function(llevel){
  require(gsw)
  tem <- subset(llevel,variable=="TEMP", value)
  psal <- subset(llevel,variable=="PSAL", value)
  depth <- llevel$depth[1] 
  
  SA <- gsw_SA_from_SP(psal,depth,llevel$lon[1],llevel$lat[1])
  CT <- gsw_CT_from_t(SA,tem,depth)
  sigma <- gsw_sigma0(SA,CT)
  rho <- gsw_rho(SA,CT,depth)  
  out <- data.frame( value = c(sigma, rho), variable = c('SIG','RHO')) #qc = rep(max(llevel$qc),2),
  return(out)
}

InSituDens <- function(llevel){
  require(gsw)
  tem <- subset(llevel,variable=="TEMP", value)
  psal <- subset(llevel,variable=="PSAL", value)
  depth <- unique(llevel$depth)
  
  SA <- gsw_SA_from_SP(psal,unique(llevel$depth),unique(llevel$lon),unique(llevel$lat))
  CT <- gsw_CT_from_t(SA,tem,depth)
  
  rho <- gsw_rho(SA,CT,depth)
  out <- data.frame( value = rho, qc = min(llevel$qc),  variable = 'RHO')
  return(out)
}


###########################################
# Functions for profiles

R20 <- function(profile){
  # Intended for Casted DF
  if (FilterForVOX(profile)){
    out <- data.frame(
    variable = "R20",
    value=min(profile$SIG[which(profile$DOXY<20)],na.rm = T)
  )
}else{
  out <- data.frame(
    variable = "Z20",
    value=NA
  )
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
    z20 <-Z20(profile)$value
    d  <- profile$depth[which(profile$depth<=z20)]
    o  <- profile$DOXY[which(profile$depth<=z20)]
    zp <- seq(1,z20)
    po <- pchip(d,o, seq(1,z20))
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
     if (is.na(out$value)|out$value>5e2|out$value<1){
      print('****')
      print(profile)
      print(min(profile$depth))
      plot(profile$DOXY,profile$depth,ylim=c(200,0) )
      lines(c(20,20),c(400,0), col='red')
      lines(pchip(d,o, seq(1,Z20)),seq(1,Z20))
      lines(po,zp,col='green')
      points(pchip(d,o, seq(1,Z20)),seq(1,Z20), pch=3)
      print(out$value)
      print('****')
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
    r  <- profile$SIG
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

FilterForVOX <- function(profile){
  # Gather here requirements for Oxygen diagnostic extraction
  flag=TRUE
  pd<-profile$depth[profile$depth<250]
  if (sum(!is.na(profile$DOXY))<10 | 
      min(profile$DOXY)>20 | 
      min(profile$depth)>15 |
      max( pd[2:length(pd)] - pd[1:(length(pd)-1)]) >40
      ){
    flag=FALSE
 
  }
  return (flag)
}

### FINISH THIS ###
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
