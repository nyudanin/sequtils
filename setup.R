#' Initial setup for RNAseq analysis
#'

#' @param project.dir Path to project directory.
#' @param nthreads Number of threads (for running featureCounts).
#' @param genome Name of reference genome assembly
#' @export
setup <- function(project.dir=getwd(),
                  nthreads=4,
                  genome="mm10"
                  ) {
  suppressPackageStartupMessages(OK <- require(Rsubread))
  if (!OK) stop("Error: Rsubread package not found")

  old.dir <- getwd()
  bam.dir <- file.path(project.dir,'data','BAM')
  bam.files <- list.files(path=bam.dir,pattern="*.bam$")

  reports.dir <- file.path(project.dir,'reports')

  dir.create(reports.dir,showWarnings=FALSE)

  setwd(bam.dir)
  capture.output(
    counts <- featureCounts(bam.files,
                          annot.inbuilt=genome,
                          useMetaFeatures=FALSE,
                          allowMultiOverlap=FALSE,
                          isPairedEnd=TRUE,
                          requireBothEndsMapped=FALSE,
                          checkFragLength=FALSE,
                          minFragLength=50,
                          maxFragLength=600,
                          nthreads=nthreads,
                          strandSpecific=0,
                          minMQS=0,
                          readExtension5=0,
                          readExtension3=0,
                          read2pos=NULL,
                          minOverlap=1,
                          countMultiMappingReads=FALSE,
                          countChimericFragments=TRUE,
                          ignoreDup=FALSE,
                          chrAliases=NULL,
                          reportReads=NULL),
    file = file.path(reports.dir,'featureCounts.report.txt')
  )

  #Set the working directory to what it was before
  setwd(old.dir)

  counts$annotation <- counts$annotation[,c('GeneID'),drop=FALSE]

  #Write out gene_counts.tsv
  write.table(x=data.frame(counts$annotation,
                           counts$counts,
                           stringsAsFactors=FALSE),
              file=file.path(project.dir,'data','gene_counts.tsv'),
              quote=FALSE,
              sep="\t",
              row.names=FALSE)
}
