OMI_figure <- function (longdf, TIMEdf, vardf, filename, addfloat = FALSE, longdf_daily=NULL ){
  library(scales)
  
  
  
trend_display <- function(df){
  m <- lm(mean ~ time, df);
  tr <- substitute( a~"+/-"~b~u,
                    list( a= format(coef(m)[2], digits = 2),
                          b = format(coef(summary(m))[,"Std. Error"][2], digits = 2), 
                          u = paste0(as.character(vardf$units),' / yr')))
                    
    as.character(as.expression(tr));                 
  }
  
  
Gplot <-
  ggplot(longdf,aes(x=as.Date(time,origin=TIMEdf$origin),y=mean))+
  geom_line(color="blue")+
  geom_ribbon(aes(ymin=mean-1.96*stde,ymax=mean+1.96*stde),alpha=0.2, fill='blue')+
  ylim(c(vardf$minsc,vardf$maxsc))+
  xlab('Year')+
  ylab(paste0('[',vardf$units,']'))+
  ggtitle(vardf$long_name)+
  theme_light()+
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")+
  geom_smooth(method = lm, color='black')+
  geom_text(x = as.Date(longdf$time[2000],origin=TIMEdf$origin),
            y = vardf$maxsc,
            label = trend_display(longdf), parse = TRUE)

if (addfloat) {
  Gplot <- Gplot + geom_point(data= melt(longdf_daily[,-which(colnames(longdf_daily) %in% c("mean","stde"))],id.vars = "time"),
                              aes(x=as.Date(time,origin=TIMEdf$origin), y = value, color=variable), 
                              alpha = 0.5)
}


print(Gplot)

pdf(width = 14, height = 10, file = filename)
print(Gplot)
dev.off()

}
