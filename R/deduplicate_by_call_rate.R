#' When a dataframe with characteristics from same inividuals (each replicate one row), for each individual, select the row with least ratio of NAs, leave each individual only one row.
#'
#' @param df a dataframe with sample "ID" as its first column
#'
#' @return a data frame with 'ID' column de-duplicated by using ID of least NA rate
#' @export
#'
#' @examples df example:
#'                ID purity ploidy
#'  TCGA-OR-A5J1-01A    0.9      2
#'  TCGA-OR-A5J2-01A   0.89    1.3
#'  TCGA-OR-A5J3-01A   0.93   1.27
#'  TCGA-OR-A5J4-01A   0.87    2.6
#'  TCGA-OR-A5J5-01A   0.93   2.79

dedup_by_callrate=function(df){
  # df=dsg_raw2

  library(magrittr)
  library(dplyr)
  library(magrittr)


  ind_info=df$ID %>% table() %>% as.data.frame() %>%
    set_colnames(c("ID","freq")) %>% mutate(ID=as.character(ID))
  uni_ind_info=ind_info %>% filter(freq==1)
  multi_ind_info=ind_info %>% filter(freq!=1)

  df_filt_raw=df %>% filter(ID %in% uni_ind_info$ID)
  for(i in multi_ind_info$ID){
    subdf=df[df$ID==i,]
    callrate_info_raw=""
    for (j in 1:NROW(subdf)) {
      CR=subdf[j,-1] %>% is.na() %>% sum()
      callrate_info_raw=rbind(callrate_info_raw,c(eachrow=j,CR=CR))
    }
    callrate_info=callrate_info_raw[-1,] %>%
      as.data.frame() %>% arrange(desc(CR)) %>%
      mutate(eachrow=as.numeric(eachrow))
    df_filt_raw=rbind(df_filt_raw,subdf[callrate_info[1,1],])
  }
  df_filt=df_filt_raw %>% remove_rownames()
  return(df_filt)
}
