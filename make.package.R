#' Create and install an R package for the project data
#'
#' @param coldata Data frame of sample data
#' @param package.name Name of R package to create
#' @export
make.package <- function(coldata,package.name){
  suppressPackageStartupMessages(OK <- require(DESeq2))
  if (!OK) stop("Error: DESeq2 package not found")
  suppressPackageStartupMessages(OK <- require(devtools))

  countdata <- read.table("data/gene_counts.tsv",
                          header=TRUE, row.names=1)
  colnames(countdata) <- rownames(coldata)
  form <- as.formula(paste("~",paste(colnames(coldata),collapse='+')))
  dds <- DESeqDataSetFromMatrix(countData=countdata,
                                colData=coldata,
                                design=form)
  dds <- dds[ rowSums(counts(dds)) > 1, ] #Get rid of any genes with 0 total counts
  dds <- DESeq(dds)
  rld <- rlogTransformation(dds)
  vsd <- varianceStabilizingTransformation(dds)

  create(package.name)
  setwd(package.name)
  use_data(coldata,countdata,dds,rld,vsd)
  setwd('..')
  install(package.name)
}
