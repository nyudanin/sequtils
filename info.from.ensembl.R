#' Returns a data frame of information (Ensembl ID, symbol and gene name) about the given genes.
#'
#' @param ensembl.vec Vector (character) of Ensembl gene identifiers
#'
#' @export
info.from.ensembl <- function(ensembl.vec){
  require(org.Mm.eg.db)
  suppressMessages(symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=ensembl.vec,keytype='ENSEMBL',column='SYMBOL',multiVals = 'list'))
  symbols.merged <- sapply(symbols,function(v) paste(v,collapse=','))
  suppressMessages(gene.names <- AnnotationDbi::mapIds(org.Hs.eg.db,keys=ensembl.vec,keytype='ENSEMBL',column='GENENAME',multiVals = 'list'))
  gene.names.merged <- sapply(gene.names,function(v) paste(v,collapse=','))
  data.frame(ENSEMBL=ensembl.vec,SYMBOL=symbols.merged,GENENAME=gene.names.merged,stringsAsFactors = FALSE)
}
