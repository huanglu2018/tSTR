
#' used to analysis the str profile, compare them to the correspinding germline profile
#'
#' @param somatic_GT rewrew
#' @param germline_GT gsdgfdgsdfg
#'
#' @return
#' @export
#'
#' @examples


G_ref="/dsk2/who/huanglu3/proj/tSTR_gtex/input/GTEX_v7_WGS_SraRunTable.txt"
T_ref="/dsk2/who/huanglu3/proj/tSTR_gtex/input/all_v7_RNAseq_SraRunTable.txt"
GTEx_WGS_manifest=fread(G_ref,data.table = F)
GTEx_RNAseq_manifest=fread(T_ref,data.table = F)

save(GTEx_WGS_manifest,file="~/proj/rpkg/tSTR/data/GTEx_WGS_manifest.RData")
save(GTEx_RNAseq_manifest,file="~/proj/rpkg/tSTR/data/GTEx_RNAseq_manifest.RData")


GT_pair_f="/dsk2/who/huanglu3/proj/tSTR_gtex/result/pair_info.rds"
GT_pair=readRDS(GT_pair_f)


str_somatic_vs_germline=function(somatic_GT,germline_GT){
  query_bed-> filted_bed



}



