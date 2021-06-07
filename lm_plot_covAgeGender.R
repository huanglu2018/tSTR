lm_plot_covAgeGender=function(Data){
  library(magrittr)
  library(ggplot2)
  
  xvar=colnames(Data)[1]
  yvar=colnames(Data)[2]
  
  
  Data=Data %>% set_colnames(c("x","y","age","gender"))
  m=lm(y~x,Data)
  coef=format(as.numeric(summary(m)$coefficients[,1][2]), digits = 3)
  r2 = format(summary(m)$r.squared, digits = 3)
  pval=format(as.numeric(summary(m)$coefficients[,4][2]), digits = 3)
  
  P=ggplot(Data,aes(x,y)) +
    geom_smooth(method = "lm",color="grey50")+
    geom_jitter(aes(color=gender,size=age),alpha=0.55,shape=20)+
    theme_bw()+labs(x = xvar, y = yvar)+
    labs(subtitle = paste0("coef = ",coef,", ","r2 = ",r2,", ","pvalue = ",pval))+
    ggsci::scale_color_nejm()
    # scale_color_gradient(low="grey80",high="black")+
    # scale_discrete_manual(values=c("#F0756C","#56BF7E"),
    #                       aesthetics = 'colour') +
    # scale_shape_manual(values=c(18,20))
    return(P)
}