#' For each individual, summarize loci infomation (DSGsum, DSGmean, DPsum, DPmean, absDSGsum, absDSGmean) by groups of loci defined by loci_list
#'
#' @param mydir dir containning .vcfs created by gangSTR
#' @param loci_list a list containing groups of loci
#'
#' @return dataframe with summarized information
#' example:
#'
#'           sample eachdis DSGsum            DSGmean DPsum           DPmean absDSGsum absDSGmean      norm_dsg norm_absdsg
#' 2 SRR2155717.vcf       1    334  0.878947368421053 16588 43.6526315789474       766  2.0157895  0.0201350374  0.04617796
#' 3 SRR2155717.vcf       2    -70 -0.171990171990172 18435 45.2948402948403       490  1.2039312 -0.0037971250  0.02657988
#' 4 SRR2155717.vcf       3     23 0.0484210526315789 21917 46.1410526315789       473  0.9957895  0.0010494137  0.02158142
#' 5 SRR2155717.vcf       4     25 0.0389408099688473 28709 44.7180685358255       491  0.7647975  0.0008708071  0.01710265
#' 6 SRR2155717.vcf       5     77  0.133913043478261 26145 45.4695652173913       629  1.0939130  0.0029451138  0.02405814
#'
#' @export
#'
#' @examples
#'
#' loci_list example as defined by exon boundary distance:
#'
#' [[1]]
#' [1] "chr1__6931916"    "chr1__10718732"   "chr1__25757550"   "chr1__31666544"   "chr1__32739601"
#'
#' [[2]]
#' [1] "chr1__36032231"   "chr1__40255942"   "chr1__40713708"   "chr1__42094312"
#'
#'
#'
get_sumabsdsg=function(mydir,loci_list){
  vcf_list=(dir(mydir) %>% grep(".vcf",.,value = T))

  res_raw=""
  count=1
  for (i in vcf_list){
    print(paste0(count,":",i))
    vcf_f=paste0(mydir,"/",i)
    if (file.info(vcf_f)$size==0) next
    DSG1=tSTR::vcf2GT(vcf_f) %>% GT2dsg() %>% mutate(loci=paste0(CHROM,"__",POS))
    DP1=tSTR::vcf2DP(vcf_f) %>% mutate(loci=paste0(CHROM,"__",POS))

    each_sum=function(idx,loci_list,DSG1,DP1,i){
      per_list=gsub("chr","",loci_list[[idx]])
      DSGsum=DSG1 %>% filter(loci%in%per_list) %>% .$DSG %>% sum()
      DSGmean=DSG1 %>% filter(loci%in%per_list) %>% .$DSG %>% mean()
      DPsum=DP1 %>% filter(loci%in%per_list) %>% .$DP %>% na.omit() %>% as.numeric() %>% sum()
      DPmean=DP1 %>% filter(loci%in%per_list) %>% .$DP %>% na.omit() %>% as.numeric() %>% mean()
      absDSGsum=DSG1 %>% filter(loci%in%per_list) %>% .$DSG %>% abs() %>% sum()
      absDSGmean=DSG1 %>% filter(loci%in%per_list) %>% .$DSG %>% abs() %>% mean()
      return(c(sample=i,eachdis=idx,
               DSGsum=DSGsum,DSGmean=DSGmean,
               DPsum=DPsum,DPmean=DPmean,
               absDSGsum=absDSGsum,absDSGmean=absDSGmean))
    }

    each_res <- foreach(idx=1:100, .combine="rbind") %do%
      {
        each_sum(idx,loci_list,DSG1,DP1,i)
      }

    res_raw=rbind(res_raw,remove_rownames(as.data.frame(each_res)))
    count=count+1
  }

  res=res_raw[-1,] %>% as.data.frame() %>%
    mutate(norm_dsg=as.numeric(DSGsum)/as.numeric(DPsum)) %>%
    mutate(eachdis=as.numeric(eachdis)) %>%
    mutate(DSGsum=as.numeric(DSGsum)) %>%
    mutate(absDSGsum=as.numeric(absDSGsum)) %>%
    mutate(absDSGmean=as.numeric(absDSGmean)) %>%
    mutate(norm_absdsg=as.numeric(absDSGsum)/as.numeric(DPsum))

  res$eachdis=as.character(res$eachdis)
  res$eachdis=factor(res$eachdis,levels = as.character(1:100))
  return(res)
}
