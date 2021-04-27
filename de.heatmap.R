#' Make a heatmap of the top differentially expressed genes
#' @param dds DESeq2 object (after running DESeq())
#' @param factor The name of the factor upon which the comparison is based (a string)
#' @param level1 The name of the first level in the comparison. Positive log fold changes correspond to higher expression in this level
#' @param level2 The name of the second level in the comparison. Negative log fold changes correspond to lower expression in this level
#' @param num.genes.up The number of upregulated genes to include in the heatmap
#' @param num.genes.down The number of downregulated genes to include in the heatmap
#' @param sample.annotations Vector of factor names to be shown in the sample annotations
#' @param cluster_cols Boolean value specifying whether to cluster columns
#' @param cluster_rows Boolean value specifying whether to cluster rows
#' @export
de.heatmap <- function(dds,factor,level1,level2,num.genes.up=50,num.genes.down=50,sample.annotations,cluster_cols=FALSE,cluster_rows=FALSE,...){
  suppressPackageStartupMessages(OK <- require(dplyr))
  if (!OK) stop("Error: dplyr package not found")
  suppressPackageStartupMessages(OK <- require(DESeq2))
  if (!OK) stop("Error: DESeq2 package not found")
  suppressPackageStartupMessages(OK <- require(pheatmap))
  if (!OK) stop("Error: pheatmap package not found")
  vsd <- varianceStabilizingTransformation(dds)
  suppressMessages(res <- ggp.results(dds,contrast=c(factor,level1,level2)))
  res %>% filter(!is.na(padj) & padj < 0.1) -> res.de
  arrange(res.de,-log2FoldChange)[1:num.genes.up,] -> top.genes.up
  arrange(res.de,log2FoldChange)[1:num.genes.up,] -> top.genes.down
  top.genes.up.entrez <- top.genes.up$ENTREZID
  top.genes.down.entrez <- top.genes.down$ENTREZID
  top.genes.up.symbol <- top.genes.up$SYMBOL
  top.genes.down.symbol <- top.genes.down$SYMBOL
  top.genes.entrez <- c(top.genes.up.entrez,top.genes.down.entrez)
  top.genes.symbol <- c(top.genes.up.symbol,top.genes.down.symbol)
  coldata <- as.data.frame(colData(dds))[,sample.annotations,drop=FALSE]
  mat <- assay(vsd)[top.genes.entrez,]
  mat <- mat - rowMeans(mat)
  rownames(mat) <- top.genes.symbol
  pheatmap(mat,
           scale='row',
           annotation_col=coldata,
           gaps_row=num.genes.up,
           show_colnames=FALSE,
           cluster_cols=cluster_cols,
           cluster_rows=cluster_rows,
           ...
    )
}
