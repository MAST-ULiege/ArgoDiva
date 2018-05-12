# ArgoVertFunctions.R

# List of function giving a scalar from an ArgoProfile

###########################################
MaxTem <- function(profile){
  tem <- subset(profile,variable=="TEMP", value)
  out <- data.frame(
    variable = "maxTem",
    value=max(tem)
  )
 return(out)
}
###########################################
MinTem <- function(profile){
  tem <- subset(profile,variable=="TEMP", value)
  out <- data.frame(
    variable = "minTem",
    value=min(tem)
  )
  return(out)
}
###########################################
MaxSal <- function(profile){
  sal <- subset(profile,variable=="PSAL", value)
  out <- data.frame(
    variable = "maxSal",
    value=max(sal)
    )
  return(out)
}
###########################################
DepthMinTem <- function(profile){
  tem <- subset(profile,variable=="TEMP")
  head(tem)
  out <- data.frame(
    variable = "depthMinTemp",
    value=mean(tem[which.min(tem$value),"depth"]) # mean in case the min tep is measured several times
  )
  return(out)
}
###########################################

DepthMinTem <- function(profile){
  tem <- subset(profile,variable=="TEMP")
  head(tem)
  out <- data.frame(
    variable = "depthMinTemp",
    value=mean(tem[which.min(tem$value),"depth"]) # mean in case the min tep is measured several times
  )
  return(out)
}
###########################################