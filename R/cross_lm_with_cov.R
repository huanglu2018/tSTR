#' cross_lm uses the variables from two data.frame and run pairwise linear regression, out put the pvalues and fdrs
#'
#' @param df1 dataframe of ncol values in nrow samples, sample ID as rowname, independent variables (x), will be lm regressed with df2 in a cross way
#' @param df2 dataframe of ncol values in nrow samples, sample ID as rowname, dependent variables (y), will be lm regressed with df1 in a cross way
#' @param df_cov dataframe of covariates
#'
#' @param fdr_method default BH, see p.adjust()
#'
#' @return var1, var2, pval, fdr, coef
#' @export
#'
#' @examples
#' df1 example:
#'                 1__100024139 1__100024569 1__100038244 1__100147494 1__100150089
#' TCGA-A4-8517-01            0           -2            0            0            0
#' TCGA-F9-A7VF-01            0           NA           NA           NA            0
#' TCGA-G7-A8LB-01           NA           NA           NA          -16            0
#' TCGA-5P-A9K3-01           -4            4            0           NA            0
#' TCGA-A4-7828-01            0            0            0           56            0
#'
#'df2 example:
#'                    aDC Adipocytes Astrocytes B-cells Basophils
#' TCGA-V4-A9EE-01 0.05877  1.375e-19 -1.621e-17 0.04805   0.06992
#' TCGA-VD-AA8N-01 0.11380  3.926e-03 -2.176e-17 0.03549   0.20740
#' TCGA-V4-A9EI-01 0.04302  3.663e-21 -1.587e-17 0.02253   0.11730
#' TCGA-VD-AA8O-01 0.08711 -2.367e-21 -2.469e-17 0.03487   0.19780
#' TCGA-WC-A888-01 0.20790 -5.244e-20 -2.434e-17 0.23040   0.18430
#'
#' cov example:
#'                  age gender
#' TCGA-2K-A9WE-01  53   male
#' TCGA-2Z-A9J1-01  71   male
#' TCGA-2Z-A9J2-01  71 female
#' TCGA-2Z-A9J3-01  67   male
#' TCGA-2Z-A9J5-01  80   male

# > df2=TIL
# > df1=dsg[,1:5]
# > df_cov=cli

cross_lm_with_cov=function(df1,df2,df_cov,fdr_method="BH"){
  options(stringsAsFactors = F)
  library(foreach)
  library(tibble)
  library(magrittr)
  library(dplyr)

  intersect_IDs<<-Reduce(intersect,list(rownames(df1),rownames(df2),rownames(df_cov)))
  print(paste0(length(intersect_IDs)," intersect IDs loaded..."))
  print(paste0(NCOL(df1)," variables loaded from dataframe 01..."))
  print(paste0(NCOL(df2)," variables loaded from dataframe 02..."))
  print(paste0(NCOL(df_cov)," variables loaded from cov..."))

  lm_var2_by_var1=function(each_var1,each_var2){
    Data=cbind(df2[intersect_IDs,each_var2,drop=F],df1[intersect_IDs,each_var1,drop=F],df_cov[intersect_IDs,]) %>%
      set_colnames(c("var2","var1",colnames(df_cov)))
    fml=as.formula(paste0("var2~var1+",paste0(colnames(df_cov),collapse = "+")))
    var2.lm=lm(fml,Data)
    pval=as.numeric(summary(var2.lm)$coefficients[,4][2])
    mycoef=as.numeric(summary(var2.lm)$coefficients[,1][2])
    return(c(var1=each_var1,var2=each_var2,pval=pval,coef=mycoef))
  }

  lm_var2_by_all_var1=function(my_var2){
    res_raw=sapply(colnames(df1),lm_var2_by_var1,each_var2=my_var2)
    res=res_raw %>% t() %>%  as.data.frame() %>%
      remove_rownames() %>%
      set_colnames(c("var1","var2","pval","coef")) %>%
      mutate(pval=as.numeric(pval)) %>%
      mutate(coef=as.numeric(coef))# %>%
    #mutate(fdr=p.adjust(pval,method = fdr_method))
    return(res)
  }

  all_res=foreach(each_var2=colnames(df2), .combine="rbind") %do%
    {
      lm_var2_by_all_var1(each_var2)
    }
  all_res=all_res %>% mutate(fdr=p.adjust(pval,method = fdr_method))
  return(all_res)
}

