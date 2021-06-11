
#' Applied to check the distance between each STR site and closest exon boundaries, if too close, it might be disturbed by alternative splicing(AS), alternative polyadenylation(APA), or alternative promoter usage(APU) events
#'
#' @param query_bed bed3 or more cols to check their overlap with references
#' @param ref_bed bed files as reference, when extracting from gtf, careful about differences between 0-index and 1-index
#'
#' @return combined df with two more columns, "inside_distance_to_ref_boundary" and "outside_distance_to_ref_boundary"; inside_distance_to_ref_boundary=0 means not inside or just at the boundary, outside_distance_to_ref_boundary=0 means not outside or just at the boundary
#' @export
#'
#' @examples
bed_overlap=function(query_bed,ref_bed){
  library(pacman)
  p_load(rtracklayer,dplyr,magrittr,tidyr)
  p_load_gh("PhanstielLab/bedtoolsr")

  bedA_raw=bt.sort(query_bed) %>% unite(.,"uniA",everything(),sep="__",remove=F)
  bedA=cbind(bedA_raw[,2:NCOL(bedA_raw)],uniA=bedA_raw[,1])

  bedB_raw=bt.sort(ref_bed) %>% unite(.,"uniB",everything(),sep="__",remove=F)
  bedB=cbind(bedB_raw[,2:NCOL(bedB_raw)],uniB=bedB_raw[,1])

  # the uniA and uniB are used for grouping

  ncA=NCOL(bedA)
  ncB=NCOL(bedB)

  res_raw=bt.closest(a = bedA, b = bedB, d=T)
  colnames(res_raw)[ncA]="uniA"
  colnames(res_raw)[ncA-3]="chrA"
  colnames(res_raw)[ncA-2]="startA"
  colnames(res_raw)[ncA-1]="endA"
  colnames(res_raw)[ncA+ncB]="uniB"
  colnames(res_raw)[ncA+1]="chrB"
  colnames(res_raw)[ncA+2]="startB"
  colnames(res_raw)[ncA+3]="endB"
  colnames(res_raw)[ncA+ncB+1]="outdistance"

  # if ref_bed is exon corrdinate, then 0 means overlap with exon, but boundary distance was not given
  res_zero=res_raw[res_raw[,ncA+ncB+1]==0,] %>%
    mutate(query_start_inside_ref=(between(startA,startB,endB))) %>%
    mutate(query_end_inside_ref=between(endA,startB,endB))
  # res_zero_inside means all the query loci was inside the corresponding ref loci
  entire_inside_exon=res_zero %>% filter(query_start_inside_ref==T & query_end_inside_ref==T) %>%
    group_by(uniA,uniB) %>%
    mutate(inside_distance_to_ref_boundary=min(abs(startA-startB),abs(startA-endB),abs(endA-startB),abs(endA-endB))) %>%
    mutate(outside_distance_to_ref_boundary=0) %>% ungroup()
  # group_by uniA means each query loci only need one distance
  # group_by uniA and uniB means each query loci may got multiple distancies to ref with ties
  only_start_inside_exon=res_zero %>%
    filter(query_start_inside_ref==T & query_end_inside_ref==F) %>%
    group_by(uniA,uniB) %>%
    mutate(inside_distance_to_ref_boundary=min(abs(startA-startB),abs(startA-endB))) %>%
    mutate(outside_distance_to_ref_boundary=min(abs(endA-startB),abs(endA-endB))) %>% ungroup()
  only_end_inside_exon=res_zero %>%
    filter(query_start_inside_ref==F & query_end_inside_ref==T) %>%
    group_by(uniA,uniB) %>%
    mutate(inside_distance_to_ref_boundary=min(abs(endA-startB),abs(endA-endB))) %>%
    mutate(outside_distance_to_ref_boundary=min(abs(startA-startB),abs(startA-endB))) %>% ungroup()

  entire_outside_exon=res_raw %>% filter(outdistance!=0) %>%
    mutate(inside_distance_to_ref_boundary=0 ) %>%
    mutate(query_start_inside_ref=F) %>%
    mutate(query_end_inside_ref=F) %>%
    group_by(uniA,uniB) %>%
    mutate(outside_distance_to_ref_boundary=min(outdistance)) %>%
    ungroup()

  res_raw=rbind(entire_inside_exon,only_start_inside_exon,only_end_inside_exon,entire_outside_exon) %>%
    select(!c(uniA,uniB,outdistance,query_start_inside_ref,query_end_inside_ref))

  res=bt.sort(res_raw)

  colnames(res)[(NCOL(res)-1):NCOL(res)]=colnames(res_raw)[(NCOL(res_raw)-1):NCOL(res_raw)]
  colnames(res)[1:(ncA-1)]=colnames(query_bed)
  colnames(res)[ncA:(ncA+ncB-2)]=colnames(ref_bed)

  return(res)
}
