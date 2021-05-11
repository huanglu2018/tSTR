


#' detect read length of fastq files, return the mean length of all the reads in given fastq files, compress status auto detected
#'
#' @param file_path fastq file, auto detect gz file
#'
#' @return mean read length
#' @export
#'
#' @examples


fq_avrage_readlen=function(file_path){

  if(!file_test("-f", file_path)) stop(simpleError('file not exist !'))
  gzcmd=paste0("gzip -dc ",file_path," | awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}'")
  cmd=paste0("awk 'NR%4==2{sum+=length($0)}END{print sum/(NR/4)}' ",file_path)
  if (endsWith(file_path,"gz")){
    len=round(as.numeric(system(gzcmd,intern = T)))
    print(paste0("read length ",len," in ", file_path))
  }else{
    len=round(as.numeric(system(cmd,intern = T)))
    print(paste0("read length ",len," in ", file_path))
  }
    return(len)
}
