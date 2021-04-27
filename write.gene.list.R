#' Write out a file with differentially expressed gene results
#'
#' @param de.genes Data frame about differentially expressed genes
#' @param filename File name for output
#' @export
write.gene.list <- function(de.genes,filename){
  write.table(de.genes,file=filename,
              sep = '\t',row.names = FALSE,quote=FALSE)
}
