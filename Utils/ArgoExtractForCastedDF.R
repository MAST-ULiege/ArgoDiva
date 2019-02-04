ArgoExtractForCastDF<- function(fulldf, flist, .parallel=FALSE){

  if (.parallel){
    cl = createCluster(8, export = c(list("flist"),flist,"FilterForVOX"), lib = list("pracma"))
  }
  
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
  }, .parallel=.parallel)#, .paropts = list(.packages=c('pracma')))
  
  if (.parallel){
    stopCluster(cl)
  }
  
  return(sdf)
}