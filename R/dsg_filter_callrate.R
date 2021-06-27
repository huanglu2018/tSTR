#' Filter dosaeg dataframe by call rate
#'
#' @param mydf  result of gather_dsg_DF, a dataframe of dosages, loci in the first column, and sampleID as colname
#' @param site_CR_threshold, anumber between 0-1, decimals
#'
#' @return
#' @export
#'
#' @examples
dsg_filter_callrate=function(mydf,site_CR_threshold){
  # use the result of gather_dsg_DF, loci as the first column
  # mydf=merge_whole_blood_RNAseq_dsg
  # site_CR_threshold=0.8
  mydf2=mydf %>% column_to_rownames("loci")
  f<-function(x) sum(is.na(x))
  site_call_rate_filt=apply(mydf2,1,f) %>%
    as.data.frame() %>%
    set_colnames("NA_num") %>%
    mutate(call_rate=1-(NA_num/NCOL(mydf2))) %>%
    filter(call_rate>=site_CR_threshold) %>%
    rownames()
  return(mydf[(mydf$loci %in%site_call_rate_filt),])
}
