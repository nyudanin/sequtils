#' @export
info.from.entrez <- function(entrez.vec){
  require(org.Hs.eg.db)
  AnnotationDbi::select(org.Hs.eg.db,
                        entrez.vec,
                        keytype="ENTREZID",
                        columns=c("ENTREZID","SYMBOL","GENENAME"))
}
