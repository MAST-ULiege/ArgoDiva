###########################################
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

###########################################
R20 <- function(profile){
  # Intended for Casted DF
  if (FilterForVOX(profile)){
    out <- data.frame(
    variable = "R20",
    value=min(profile$RHO[which(profile$DOXY<20)],na.rm = T)
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
  if (profile$Platform=="6900804 "|out$value<20 |out$value>300){
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
  if (FilterForCCC(profile)){
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
      variable = "CCC",
      value=NA
    )
  }
  return(out)
}





FilterForVOX <- function(profile){
  # Gather here requirements for Oxygen diagnostic extraction
  flag=TRUE
  if (sum(!is.na(profile$DOXY))<10 | min(profile$DOXY)>20 | min(profile$depth)>15){
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
