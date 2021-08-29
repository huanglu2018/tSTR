#' multivariate linear regression using gender and age as covariates
#'
#' @param Data dataframe with FOUR columns, x,y,age,gender. formula: y~x+age+gender, exclude gender if gender are the same
#'
#' @return a pic with coef, adj.r.squared, and pvalue
#' @export
#'
#' @examples
#'


lm_plot_covAgeGender=function(myData){
  library(magrittr)
  library(ggplot2)

  xvar=colnames(myData)[1]
  yvar=colnames(myData)[2]

  myData=myData %>% set_colnames(c("x","y","age","gender"))
  if(length(unique(myData$gender))<2) {
    print("all sample with same gender, exclude gender from regression")
    m=lm(y~x+age,myData)
  }else{
    m=lm(y~x+age+gender,myData)
  }
  coef=format(as.numeric(summary(m)$coefficients[,1][2]), digits = 3)
  r2 = format(summary(m)$adj.r.squared, digits = 3)
  pval=format(as.numeric(summary(m)$coefficients[,4][2]), digits = 3)

  P=ggplot(myData,aes(x,y)) +
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
