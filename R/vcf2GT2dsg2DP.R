
#' extract infomation of STR from gangSTR vcfs
#'
#' @param vcf_file the out file of gangSTR
#'
#' @return GT, a dataframe extracted from eachvcf
#' @export NA
#'
#' @examples
#' s=vcf2GT(vcf_file) %>% GT2maxlen()
#' s=vcf2GT(vcf_file) %>% GT2dsg()
#' s=vcf2DP(vcf_file)
#'
#' GT
#' ChromKey,POS,Indiv,gt_GT,gt_GT_alleles,REF
#' 1,14070,SRR656421Aligned.sortedByCoord.out,.,,cctccctccctc
#' 1,16620,SRR656421Aligned.sortedByCoord.out,0/0,gctgctgctgct/gctgctgctgct,gctgctgctgct
#' 1,22812,SRR656421Aligned.sortedByCoord.out,0/0,aggaaaggaa/aggaaaggaa,aggaaaggaa
#' 1,26454,SRR656421Aligned.sortedByCoord.out,.,,gtgtgtgtgtgt
#'
#' dsg:
#' ChromKey    POS SRR600493Aligned.sortedByCoord.out
#'     1  16620                                  0
#'     1  22812                                 10
#'     1 138589                                -12
#'
#'
#' maxlen:
#'     1  16620                                 12
#'     1  22812                                 15
#'     1 138589                                 12
#'
#' DP:
#'  CHROM   POS   DP
#'  chr1 14070 <NA>
#'  chr1 16620    6
#'  chr1 22812    2
#'
#'

vcf2GT=function(vcf_file){
  pacman::p_load(magrittr,data.table,dplyr)
  if(!file_test("-f", vcf_file)) stop(simpleError(paste0(vcf_file,' not exist !')))
  s=fread(vcf_file,data.table = F,skip="#CHROM") %>%
    .[,c(1,2,4,5,10)]
  s$Indiv=colnames(s)[5]
  s0=s %>%
    set_colnames(c("ChromKey","POS","REF","ALT","gt_GT","Indiv")) %>%
    mutate(ChromKey=sub("^chr","",ChromKey)) %>%
    mutate(ALT=ifelse(ALT==".",REF,ALT)) %>%
    mutate(gt_GT=substr(gt_GT,1,3)) %>%
    mutate(gt0_0=paste0(REF,"/",REF)) %>%
    mutate(gt0_1=paste0(REF,"/",ALT)) %>%
    mutate(gt1_0=paste0(ALT,"/",REF)) %>%
    mutate(gt1_2=sub(",","/",ALT)) %>%
    mutate(gt1_1=paste0(ALT,"/",ALT))
  s0[s0$gt_GT=="0/0","gt_GT_alleles"]=s0[s0$gt_GT=="0/0","gt0_0"]
  s0[s0$gt_GT=="0/1","gt_GT_alleles"]=s0[s0$gt_GT=="0/1","gt0_1"]
  s0[s0$gt_GT=="1/0","gt_GT_alleles"]=s0[s0$gt_GT=="1/0","gt1_0"]
  s0[s0$gt_GT=="1/2","gt_GT_alleles"]=s0[s0$gt_GT=="1/2","gt1_2"]
  s0[s0$gt_GT=="1/1","gt_GT_alleles"]=s0[s0$gt_GT=="1/1","gt1_1"]

  return(s0[,c("ChromKey","POS","Indiv","gt_GT","gt_GT_alleles","REF")])
}


GT2dsg=function(GT){
  pacman::p_load(tidyr,plyr,dplyr)
  s1=GT %>% filter(gt_GT%in%c("0/0","0/1","1/0","1/2","1/1")) %>%
    separate(.,5,c("allel1","allel2"),"/") %>%
    mutate(len1=nchar(allel1)) %>%
    mutate(len2=nchar(allel2)) %>%
    mutate(len3=nchar(REF)) %>%
    mutate(dsg=as.numeric(len1)+as.numeric(len2)-(2*as.numeric(len3))) %>%
    .[,c(1,2,11)] %>% set_colnames(c("CHROM","POS","DSG"))
  return(s1)
}


GT2maxlen=function(GT){
  pacman::p_load(tidyr,plyr,dplyr)
  s1=GT %>% filter(gt_GT%in%c("0/0","0/1","1/0","1/2","1/1")) %>%
    separate(.,5,c("allel1","allel2"),"/") %>%
    mutate(len1=nchar(allel1)) %>%
    mutate(len2=nchar(allel2)) %>%
    mutate(len3=nchar(REF))
  s2=data.frame(s1, maxlen = apply(s1[,c("len1","len2","len3")],1,max)) %>%
    .[,c(1,2,11)] %>% set_colnames(c("CHROM","POS","MAXLEN"))
  return(s2)
}


vcf2DP=function(vcf_file){
  pacman::p_load(magrittr,data.table,dplyr)
  if(!file_test("-f", vcf_file)) stop(simpleError(paste0(vcf_file,' not exist !')))
  s=fread(vcf_file,data.table = F,skip="#CHROM",select = c(1,2,10))
  DP=s[,3] %>% tstrsplit(.,":") %>% as.data.frame() %>% .[,2,drop=F]
  s2=data.frame(s[,1:2],DP) %>% set_colnames(c("CHROM","POS","DP")) %>% mutate(CHROM=sub("^chr","",CHROM))
  return(s2)
}


#vcf_file="/public5/huanglu/gtex2tSTR/out/skin-sun-exposed/result/gangstr/SRR600493.vcf"
# DP_example=vcf2DP(vcf_file)
# GT_example=vcf2GT(vcf_file)
# DSG_example=vcf_file %>% vcf2GT() %>% GT2dsg()
# MAXLEN_example=vcf_file %>% vcf2GT() %>% GT2maxlen()
