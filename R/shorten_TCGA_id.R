#' used to conveniently cut TCGA sample IDs
#'
#' @param id_list individual list, e.g. TCGA-AAAA-BBBB
#' @param left_number how many parts (separated by '-' to keep, should not be longer than its original length )
#'
#' @return a shortened list
#' @export
#'
#' @examples
#'


shorten_TCGA_id=function(id_list,left_number){
  library(data.table)
  library(tidyr)
  library(magrittr)
  library(dplyr)

  s=id_list %>% as.data.frame() %>%
    set_colnames(c("all")) %>%
    separate(.,1,paste0("V",1:7),"-",fill="right") %>%
    unite(x, c(paste0("V",1:as.numeric(left_number))), sep = "-", remove = FALSE) %>%
     .$x
  return(s)
}
