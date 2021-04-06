#' Title
#'
#' @param Data dataframe with two columns, x column and y column
#'
#' @return a pic with coef, r2, and pvalue
#' @export
#'
#' @examples
#'

lm_plot=function(Data){

  library(magrittr)
  library(ggplot2)

  xvar=colnames(Data)[1]
  yvar=colnames(Data)[2]

  Data=Data %>% set_colnames(c("x","y"))
  m=lm(y~x,Data)
  coef=format(as.numeric(summary(m)$coefficients[,1][2]), digits = 3)
  r2 = format(summary(m)$r.squared, digits = 3)
  pval=format(as.numeric(summary(m)$coefficients[,4][2]), digits = 3)

  P=ggplot(Data,aes(x,y)) +
    geom_jitter()+
    geom_smooth(method = "lm")+
    theme_bw()+labs(x = xvar, y = yvar)+
    labs(subtitle = paste0("coef = ",coef,", ","r2 = ",r2,", ","pvalue = ",pval))
  return(P)
}

