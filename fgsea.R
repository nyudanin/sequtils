#' Run Fast GSEA (using fgsea package by A. A. Sergushichev)
#'
#' @param dds DESeq2 object to use
#' @param attr Name (character) of sample attribute to use in the comparison, e.g. 'Genotype'
#' @param level1 Name (character) of first level in comparison, e.g. 'Knockout'. Positive values of the test statistic will correspond to higher expression in this level.
#' @param level2 Name (character) of second level in comparison.
#' @param genesets A list vectors of gene Entrez IDs (type: character) to be used in fgsea. By default, all GO terms are used.
#' @param minsize The minimum number of genes in a gene set in order for it to be included in the fgsea analysis. Default: 5
#' @param maxsize The maximum number of genes in a gene set in order for it to be included in the fgsea analysis. Default: 5000
#' @param nperm The number of permutations used to calculate p values. Default: 10,000
#' @export
run.fgsea <- function(dds,attr,level1,level2,genesets='GO',minsize=5,maxsize=5000,nperm=10000){
  suppressPackageStartupMessages(OK <- require(fgsea))
  if (!OK) stop("Error: package 'fgsea' not found")
  if (genesets=='GO'){
    use.go <- TRUE
    suppressPackageStartupMessages(OK <- require(ggp.genesets))
    if (!OK) stop("Error: package 'ggp.genesets' not found")
    genesets <- all.go
  } else{
    use.go <- FALSE
  }
  suppressMessages(res <- results(dds,contrast=c(attr,level1,level2)))
  ord <- order(-res$stat)
  res.ranked <- res[ord,] #Results are ranked by the test statistic
  ranks <- res.ranked$stat
  names(ranks) <- res.ranked$ENTREZID

  #Run fgsea
  fgseaRes <- fgsea(pathways=genesets,stats=ranks,minSize=minsize,maxSize=maxsize,nperm=nperm)
  if (use.go){
    #Include GO pathway name as well as GOID, and rename the 'pathway' column to 'GOID'
    colnames(fgseaRes)[1] <- 'GOID'
    fgseaRes <- data.frame(name=all.go.names[fgseaRes$GOID],fgseaRes)
  }
  fgseaRes
}
