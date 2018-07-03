ArgoComplete<- function(fulldf, flist){
  
  sdf <- ddply(fulldf, ~ aprofile+Platform+alevel, function(llevel){
    flistoutputs <- ldply(flistl,function(f){
      f<- match.fun(f)
      outf<-f(llevel)
      }
    )
    
#    bi<-rep(llevel[1,],nrow(flistoutputs))
    bi<-llevel[1,]
    bi$value<-flistoutputs$value
    bi$variable<-flistoutputs$variable
    
    return( rbind(llevel, bi) )
  })
  
  
  return(sdf)
}