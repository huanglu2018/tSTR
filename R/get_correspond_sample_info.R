#' Function used to map GTEx v7 SRR id(s) from RNAseq to its tissue part, WGS SRR id and individual id
#'
#'
#' @param inSRR one or a list of GTEx RNAseq SRR id
#' @param outdf a dataframe idicating the mapping result
#'
#' @return a dataframe containing columns:individual,RNAseq_run,RNAseq_sample,RNAseq_tissue,WGS_run,WGS_sample,WGS_tissue
#' @export
#'
#' @examples
#'
#' inSRR:
#' c("SRR598300")
#'
#' out_df:
#'   individual   RNAseq_run   RNAseq_sample   RNAseq_tissue   WGS_run   WGS_sample   WGS_tissue
#'   GTEX-PX3G   SRR598300   GTEX-PX3G-0006-SM-2I3E4   Whole Blood   SRR3478911   GTEX-PX3G-0004-SM-5JK53   Whole Blood
#'
get_correspond_sample_info=function(inSRR){
  pacman::p_load(magrittr,data.table,dplyr,tibble)
  # inSRR=c("SRR598300")
  data("GTEx_RNAseq_manifest")
  data("GTEx_WGS_manifest")
  T_info=GTEx_RNAseq_manifest[,c("Run","Sample_Name","body_site","submitted_subject_id")] %>%
    set_colnames(c("RNAseq_run","RNAseq_sample","RNAseq_tissue","RNAseq_individual")) %>%
    filter(RNAseq_run %in% inSRR)
  G_info=GTEx_WGS_manifest[,c("Run","Sample_Name","body_site","submitted_subject_id")] %>%
    set_colnames(c("WGS_run","WGS_sample","WGS_tissue","WGS_individual")) %>%
    filter(WGS_individual %in% T_info$RNAseq_individual)
  outdf=merge(T_info,G_info,by.x="RNAseq_individual",by.y="WGS_individual")
  colnames(outdf)[1]="individual"
  rownames(outdf)=outdf$RNAseq_run
  outdf2=outdf[inSRR,] %>% remove_rownames()
  return(outdf2)
}
