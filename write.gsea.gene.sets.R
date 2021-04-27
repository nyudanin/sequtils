#' Create .gmt file of mouse GO terms
#'
#' @param ontologies Vector of GO ontologies to be included: a subset of c('BP','MF','CC').
#' @param filename Output filename, typically ending in ".gmt".
#' @export
write.gsea.gene.sets.go <- function(ontologies='BP',filename)
{
  suppressPackageStartupMessages(OK <- require(ggp.genesets))
  if (!OK) stop("Error: package 'ggp.genesets' not found")
  suppressPackageStartupMessages(OK <- require(GO.db))
  if (!OK) stop("Error: package 'GO.db' not found")

  go.strs <- sapply(all.go,function(x) paste(x,collapse='\t'))
  go.data <- data.frame(
    row.names=names(all.go),
    term=suppressMessages(AnnotationDbi::mapIds(GO.db,keys = names(all.go),keytype='GOID',column='TERM')),
    ontology=suppressMessages(AnnotationDbi::mapIds(GO.db,keys = names(all.go),keytype='GOID',column='ONTOLOGY')),
    strs=go.strs
  )
  go.data.filtered <- go.data[go.data$ontology %in% ontologies,c('term','strs')]

  go.data.filtered.relabelled <- go.data.filtered
  rownames(go.data.filtered.relabelled) <- paste(rownames(go.data.filtered),go.data.filtered$term,sep=' ')
  go.data.filtered.relabelled$term <- rownames(go.data.filtered)

  go.data.filtered.relabelled.nospace <- go.data.filtered.relabelled
  rownames(go.data.filtered.relabelled.nospace) <- sapply(rownames(go.data.filtered.relabelled.nospace),
                                                    function(x) gsub(' ','_',x))

  write.table(go.data.filtered.relabelled.nospace,
              file=filename,sep='\t',quote=FALSE,col.names=FALSE)
}
