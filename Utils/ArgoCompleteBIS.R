ArgoCompleteBIS<- function(fulldf, flistl, .parallel=FALSE){
  
  if (.parallel){
  cl = createCluster(8, export = c(list("flistl"),flistl), lib = list("gsw"))
  }
  
  sdf <- ddply(fulldf, ~ juld+Platform+depth, function(llevel){
  
    flistoutputs <- ldply(flistl,function(f){
      f<- match.fun(f)
      outf<-f(llevel)
      }
    )
    bi<-llevel[1,]
    bi <- bi[rep(1,length(flistoutputs$value)),]
    bi$value<-flistoutputs$value
    bi$variable<-flistoutputs$variable
    
    return( rbind(llevel, bi) )
  }, .parallel=.parallel)#, .paropts = list(.packages=c('gsw'),.export="flistl"))
  
  if (.parallel){
  stopCluster(cl)
  }
  return(sdf)
}
