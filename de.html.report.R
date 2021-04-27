#' HTML report for a differential expression comparison
#'
#' @param dds DESeq object
#' @param output.dir Output directory name to be created for the HTML report
#' @param attr The attribute to be compared (i.e. 'genotype' or 'treatment')
#' @param level1 The value of the comparison attribute in group 1, which will be on the y axis (i.e. 'Knockout' or 'Infected')
#' @param level2 The value of the comparison attribute in group 2, which will be on the x axis (i.e. 'Wildtype' or 'Healthy')
#' @param alpha Significance cutoff for highlighting genes
#' @param annotation.db The annotation database (i.e. \code{org.Mm.eg.db}) used to look up information about features (genes or transcripts)
#' @param identifier The name of the column of \code{annotation.db} that corresponds to the row names of \code{dds} (feature names)
#' @param gene.info.columns Vector of column names from \code{annotation.db} to include in the output (for example, \code{'SYMBOL'}, \code{'GENENAME'}...)
#' @export
de.html.report <- function(dds,output.dir,attr,level1,level2,alpha=0.1,annotation.db=NULL,identifier=NULL,gene.info.columns=NULL){
  args.list <- as.list(environment(), all=TRUE)

  # First, get the annotation database and the identifier
  annotation.info <- annotation.db.helper(args.list,obj=dds)

  res.df <- ggp.results(dds,contrast=c(attr,level1,level2),annotation.db=annotation.db,identifier=identifier,gene.info.columns=gene.info.columns)

  html.rep <- ReportingTools::HTMLReport(
    shortName=paste0(c(level1,'_vs_',level2),collapse=''),
    title=paste0(c('Differential expression results: ',level1, ' vs. ',level2),collapse=''),
    reportDirectory = output.dir
  )
  ReportingTools::publish(paste0("Positive values of logFC indicate higher expression in ",level1,"."),html.rep)
  ReportingTools::publish(paste0("Analysis package: DESeq2 version ",metadata(dds)$version),html.rep)
  ReportingTools::publish(paste0("Statistical test: ",attr(dds,'test')),html.rep)
  ReportingTools::publish("p value adjustment: False Discover Rate (FDR)",html.rep)
  ReportingTools::publish(paste0("Criterion for being called differentially expressed: adjusted p < ",alpha),html.rep)
  mdata <- attr(res.df,'metadata')
  ReportingTools::publish(paste0("Independent filtering threshold: ",as.character(mdata$filterThreshold)),html.rep)
  res.diff.exp <- res.df[!is.na(res.df$padj) & res.df$padj < 0.1,]

  if(sum(!is.na(res.df$padj) & res.df$padj < 0.1) > 0){
    ReportingTools::publish(res.diff.exp,html.rep)
  }
  else{
    ReportingTools::publish("<B><H2>No differentially expressed genes.</H2></B>",html.rep)
  }
  ReportingTools::finish(html.rep)
}
