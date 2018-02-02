# ArgoVertFunctions.R

# List of function giving a scalar from an ArgoProfile

MaxTem <- function(profile){
  tem <- subset(profile,variable=="TEMP", value)
  out <- data.frame(
    variable = "maxTem",
    value=max(tem)
  )
 return(out)
}

MaxSal <- function(profile){
  sal <- subset(profile,variable=="PSAL", value)
  out <- data.frame(
    variable = "maxSal",
    value=max(sal)
    )
  return(out)
}
