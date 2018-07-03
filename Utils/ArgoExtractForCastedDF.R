require("plyr")

ArgoExtractForCastDF<- function(fulldf, flist){

  sdf <- ddply(fulldf, ~ juld+Platform, function(profile){
    flistoutputs <- ldply(flist,function(f){
      f<- match.fun(f)
      f(profile)
    }
    )
    
    onelinedf <- data.frame(
      qc    = max(profile$qc), # TODO  several variables , etc ... 
      lon   = mean(profile$lon),
      lat   = mean(profile$lat) , 
      day   = profile$day[1], 
      month = profile$month[1], 
      year  = profile$year[1],
      juld  = profile$juld[1])
    
    return( cbind(flistoutputs, onelinedf) )
  })
  
  
  return(sdf)
}