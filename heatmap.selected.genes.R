#' Make a heatmap from selected genes
#' @param dds.obj The DESeq2 object
#' @param genes Vector of gene identifiers
#' @param sample.annotations Vector of factor names to be shown in the sample annotations
#' @param identifier The type of gene identifier ('ENTREZID','SYMBOL'). Default is 'ENTREZID'.
#' @export
heatmap.selected.genes <- function(dds.obj,
                                       genes,
                                       sample.annotations,
                                       identifier='ENTREZID',
                                       ...
                                      ){
  suppressPackageStartupMessages(OK <- require(DESeq2))
  if (!OK) stop("Error: DESeq2 package not found")
  suppressPackageStartupMessages(OK <- require(pheatmap))
  if (!OK) stop("Error: pheatmap package not found")
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  if (!OK) stop("Error: org.Hs.eg.db package not found")

  switch(identifier,
         ENTREZID={
           entrez <- genes
           symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                           keys=entrez,
                                           keytype='ENTREZID',
                                           column="SYMBOL",
                                           multiVals='first'
           )
         },
         SYMBOL={
          symbol <- genes
          entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                          keys=genes,
                                          keytype='SYMBOL',
                                          column="ENTREZID",
                                          multiVals='first'
          )
         },
         { #Default
          stop("Error: invalid identifier type: ", identifier)
         }
        )
  vsd <- varianceStabilizingTransformation(dds.obj)

  common <- intersect(entrez,rownames(assay(vsd)))
  not.included <- setdiff(entrez,rownames(assay(vsd)))

  common.symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                  keys=common,
                                  keytype='ENTREZID',
                                  column="SYMBOL",
                                  multiVals='first'
  )

  if(length(not.included) > 0){
    message("The following genes were not found in the expression data: ", paste(not.included,collapse=', '))
  }

  mat <- assay(vsd)[common,,drop=FALSE]
  mat <- mat - rowMeans(mat)
  rownames(mat) <- common.symbol

  coldata <- as.data.frame(colData(dds))[,sample.annotations,drop=FALSE]

  pheatmap(mat,
           scale='row',
           annotation_col=coldata,
           show_colnames=FALSE,
           cluster_cols=FALSE,
           cluster_rows=FALSE,
           ...
  )
}
