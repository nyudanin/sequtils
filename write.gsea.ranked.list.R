#' Create a ranked gene list for use with GSEA's GseaPreranked mode
#'
#' @param dds DESeq2 object to be used for generating the ranked list.
#' @param attr The attribute to be compared (i.e. 'genotype' or 'treatment')
#' @param level1 The value of the comparison attribute in group 1 (i.e. 'Knockout' or 'Infected')
#' @param level2 The value of the comparison attribute in group 1 (i.e. 'Wildtype' or 'Healthy')
#' @param ranking.metric The name of the metric by which to rank genes. Must be a column of DESeq2's `results` output. Default: "stat".
#' @param identifier The type of gene identifier: "ENTREZID" or "SYMBOL". Default: "ENTREZID"
#' @param filename The name of the output file (typically ending with ".rnk")
#' @export
write.gsea.ranked.list <- function(dds,attr,level1,level2,ranking.metric='stat',identifier='ENTREZID',filename){
  suppressPackageStartupMessages(OK <- require(dplyr))
  if (!OK) stop("Error: package 'dplyr' not found")

  res <- ggp.results(dds,contrast=c(attr,level1,level2))
  res %>% arrange_(ranking.metric) -> res.ranked

  write.table(res.ranked[,c(identifier,ranking.metric)],
              file=filename,
              quote=F,
              sep="\t",
              row.names=F)
}
