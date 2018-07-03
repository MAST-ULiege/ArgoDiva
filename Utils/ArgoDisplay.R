 ArgoDisplay<- function(fdf,selectCriterium,flist){

   rmarkdown::render("ArgoDiagReport.Rmd", params = list(
       fdf              = fdf  ,
       selectCriterium  = selectCriterium , 
       flist            = flist
     ))
     
   
 
   extdir='./Reports/'
  dir.create(extdir)
   file.copy(from = paste0(getwd(),"/ArgoDiagReport.pdf"),
               to = paste0(extdir,"ArgoDiagReport_",Sys.Date(),".pdf"), overwrite = T)
   
   
    
 }
