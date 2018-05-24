#' Extract particular fields from the attribute column in the gff
#' 
#' @param gff a vector where values are the attributes in GFF2 format
#' @param attribute the name of the attribute to extract
getGff2Attribute <- function(gff, attribute){
  str_extract(gff, paste0(attribute, " .*")) %>% 
    str_remove(" ;.*") %>% 
    str_remove(paste0(attribute, " "))
}
