

vcf_extract_GT=function(vcf_file){
  library(magrittr)
  library(data.table)
  library(dplyr)
  s=fread(vcf_file,data.table = F,skip="#CHROM") %>%
    .[,c(1,2,4,5,10)]
  s$Indiv=gsub("Aligned.sortedByCoord.out$","",colnames(s)[5])
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

  s0[,c("ChromKey","POS","Indiv","gt_GT","gt_GT_alleles","REF")]
}

get_dsg=function(GT,read_len){
  GT %>% separate(.,5,c("allel1","allel2"),"/") %>%
    mutate(len1=nchar(allel1)) %>%
    mutate(len2=nchar(allel2)) %>%
    mutate(len3=nchar(REF)) %>%
    ddply(.,.(ChromKey,POS),summarise,
          maxlen=max(len1,len2,len3),Indiv=Indiv,gt_GT=gt_GT,
          allel1=allel1,allel2=allel2,REF=REF,
          len1=len1,len2=len2,len3=len3) %>%
    mutate(maxlen=ifelse(is.na(maxlen),len3,maxlen)) %>%
    mutate(dsg=ifelse(maxlen<=read_len,
                      as.numeric(len1)+as.numeric(len2)-(2*as.numeric(len3)),NA)) %>%
    .[,c(1,2,12)]
}



ganstrvcf2dosage=function(vcf_file,res_GT_file,res_dsg_file){
  each_GT=try(vcf_extract_GT(vcf_file))
  if ('try-error' %in% class(each_GT)) {
    print(paste0("GT-extract error in ",i))
    next
  }else{
    fwrite(each_GT,file = each_GT_f,compress = "gzip")
    each_read_len=len_info[match(ID,len_info[,1]),2]
    if(file.exists(each_dsg_f)) next
    each_dsg=try(get_dsg(each_GT,each_read_len))
    if ('try-error' %in% class(each_dsg)) {
      print(paste0("dosage-extract error in ",i))
      next
    }else{
      fwrite(each_dsg,file = each_dsg_f,compress = "gzip")
    }
  }
}




if(!dir.exists(res_dsg_dir)) dir.create(res_dsg_dir,recursive = T)
if(!dir.exists(re_GT_dir)) dir.create(re_GT_dir,recursive = T)
len_info=fread(len_info_file,data.table = F)

for (i in c(grep(".vcf$",dir(vcf_dir),value=T),grep(".vcf.gz$",dir(vcf_dir),value=T))){
  ID=gsub(".vcf.gz$","",i) %>% gsub("
                                    $","",.)
  each_dsg_f=paste0(res_dsg_dir,"/",ID,".dsg.gz")
  each_GT_f=paste0(re_GT_dir,"/",ID,".gt.gz")
  if(file.exists(each_GT_f)) next
  print(paste0("====== working on ",i," ... ======"))

  each_GT=try(vcf_extract_GT(paste0(vcf_dir,"/",i)))
  if ('try-error' %in% class(each_GT)) {
    print(paste0("GT-extract error in ",i))
    next
  }else{
    fwrite(each_GT,file = each_GT_f,compress = "gzip")
    each_read_len=len_info[match(ID,len_info[,1]),2]
    if(file.exists(each_dsg_f)) next
    each_dsg=try(get_dsg(each_GT,each_read_len))
    if ('try-error' %in% class(each_dsg)) {
      print(paste0("dosage-extract error in ",i))
      next
    }else{
      fwrite(each_dsg,file = each_dsg_f,compress = "gzip")
    }
  }
}



# 计算dosage的条件，两个等位位点的长度减去2*REF的长度，如果这三个的长度有一处长于测序的read_len,则dsg也为NA
# 有的vcf里面虽然REF和ALT都非空，但是最后一列的GT为“.”，这种也没法用





