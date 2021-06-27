
#' used to analysis the str profile, compare them to the correspinding germline profile
#'
#' @param RNAseq_GT result of tSTR::vcf2GT() and tSTR::get_correspond_sample_info()
#' @param WGS_GT result of tSTR::vcf2GT() and tSTR::get_correspond_sample_info()
#'
#' @return a dataframe include columns:loci T_gt_GT,TG1,TG2,T_REF,G_gt_GT,GG1,GG2,G_REF,dsg_modi,dsg_mode
#' @export
#'
#' @examples
#'
#' result dataframe
#'
#'       loci    T_gt_GT  TG1    TG2      T_REF         G_gt_GT       GG1           GG2           G_REF            dsg_modi   dsg_mode
#'  1__10034928     1/1   gatg   gatg gatggatggatggatg     0/0 gatggatggatggatg gatggatggatggatg gatggatggatggatg       24 mutate_all
#'  1__10035556     1/1  aaaac  aaaac  aaaacaaaacaaaac     0/0  aaaacaaaacaaaac  aaaacaaaacaaaac  aaaacaaaacaaaac       20 mutate_all
#'  1__100405174     1/1 ttgttg ttgttg  ttgttgttgttgttg     0/0  ttgttgttgttgttg  ttgttgttgttgttg  ttgttgttgttgttg       18 mutate_all
#'



#
# dirG="/dsk2/who/huanglu3/data/vcf/bed_131_vcf_image_result"
# dirT="/public4/huanglu/gtex2tSTR/out/whole-blood/result/gangstr"
#
# inSRR=c("SRR1488261")
# T_path=paste0(dirT,"/",inSRR,".vcf")
# G_path=paste0(dirG,"/",get_correspond_sample_info(inSRR)$WGS_run,".vcf")
# file.exists(T_path)
# file.exists(G_path)
#
# T_path
# G_path
#
# RNAseq_GT=tSTR::vcf2GT(T_path)
# WGS_GT=tSTR::vcf2GT(G_path)

str_somatic_vs_germline=function(RNAseq_GT,WGS_GT){
  pacman::p_load(magrittr,data.table,dplyr,tibble,foreach)
  T1=RNAseq_GT %>%
    na.omit() %>%
    mutate(loci=paste0(ChromKey,"__",POS)) %>%
    .[,c(7,4,5,6)] %>%
    set_colnames(c("loci","T_gt_GT","T_gt_GT_alleles","T_REF"))

  G1=WGS_GT %>%
    na.omit() %>%
    mutate(loci=paste0(ChromKey,"__",POS)) %>%
    .[,c(7,4,5,6)] %>%
    set_colnames(c("loci","G_gt_GT","G_gt_GT_alleles","G_REF"))

  M_dsg=merge(T1,G1,by="loci") %>%
    separate(.,3,c("TG1","TG2"),"/") %>%
    separate(.,7,c("GG1","GG2"),"/")

  each_row_compare=function(i, M_dsg){
    GG1=M_dsg[i,"GG1"]
    GG2=M_dsg[i,"GG2"]
    TG1=M_dsg[i,"TG1"]
    TG2=M_dsg[i,"TG2"]
    loci=M_dsg[i,"loci"]

    if ((GG1==TG1 & GG2==TG2)|(GG1==TG2 & GG2==TG1)){
      dsg_modi=0
      dsg_mode="same"
    }else if(GG1==TG1 | GG2==TG2 | GG1==TG2 | GG2==TG1){
      dsg_modi=nchar(GG1)+nchar(GG2)-nchar(TG1)-nchar(TG2)
      dsg_mode="mutate1"
    }else{
      dsg_modi=nchar(GG1)+nchar(GG2)-nchar(TG1)-nchar(TG2)
      dsg_mode="mutate_all"
    }
    return(c(loci=loci,dsg_modi=dsg_modi,dsg_mode=dsg_mode))
  }

  modi_info <- foreach(i=seq(NROW(M_dsg)), .combine="rbind") %do%
    {
      each_row_compare(i,M_dsg)
    }

  M_dsg2=merge(M_dsg,as.data.frame(modi_info),by="loci")
  return(M_dsg2)
}



