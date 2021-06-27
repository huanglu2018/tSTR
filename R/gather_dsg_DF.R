#' Gather dosage information to a dataframe from a directory containing gangSTR vcfs (wrapper of vcf2GT and GT2dsg)
#'
#' @param mydir the directory of vcf files created by gangSTR
#'
#' @return a dataframe of dosages, loci in the first column, and sampleID as colname
#' @export
#'
#' @examples
gather_dsg_DF=function(mydir){
  wb_dir=mydir
  vcf_list=(dir(wb_dir) %>% grep(".vcf",.,value = T))

  merge_dsg_f=paste0(wb_dir,"/",vcf_list[1])
  merge_dsg=tSTR::vcf2GT(merge_dsg_f) %>%
    tSTR::GT2dsg() %>%
    mutate(loci=paste0(CHROM,"__",POS)) %>%
    mutate(dummy_dsg=DSG) %>%
    .[,4:5]

  count=1
  for (i in vcf_list){
    print(paste0(count,":",i))
    vcf_f=paste0(wb_dir,"/",i)
    if (file.info(vcf_f)$size==0) next
    dsg_i=tSTR::vcf2GT(vcf_f) %>%
      tSTR::GT2dsg() %>%
      mutate(loci=paste0(CHROM,"__",POS)) %>%
      .[,3:4]
    merge_dsg_raw=merge(merge_dsg,dsg_i,by="loci",all = T)
    merge_dsg=plyr::rename(merge_dsg_raw,c("DSG"=gsub(".vcf","",i)))
    count=count+1
  }
  return(merge_dsg[,-2])
}
