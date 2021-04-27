#' Volcano plot
#'
#' Make a volcano plot from the results dataframe \code{res.df}, which should be the output of `ggp.results`
#' Significantly differentially expressed genes (using adjusted p value threshold \code{alpha}) are highlighted in red.
#' Genes will be labelled as specified by \code{labelled.genes}, which must be identifiers of type \code{labelled.genes.keytype}.
#' As in many \code{ggp.rnaseq} functions, \code{annotation.db} can be supplied via the \code{'annotation.db'} attribute of the
#' object \code{dds}, and will default to \code{org.Mm.eg.db}. Once the genes
#' \code{labelled.genes} have been selected, they are labelled in the plot using the column \code{display.identifier} of \code{annotation.db}.
#' If \code{display.identifier} is not supplied, \code{labelled.genes.keytype} is used.
#'
#' @param res.df DE results data frame
#' @param alpha Significance cutoff for highlighting genes
#' @param labelled.genes List of genes/features to be labelled
#' @param labelled.genes.keytype The gene identifier type (e.g. \code{'ENTREZID'}) of \code{labelled.genes}
#' @param annotation.db The annotation database (i.e. \code{org.Mm.eg.db}) used to look up information about features (genes or transcripts)
#' @param identifier The name of the column of \code{annotation.db} that corresponds to the row names of \code{res.df} (feature names)
#' @param gene.info.columns Vector of column names from \code{annotation.db} to include in the output (for example, \code{'SYMBOL'}, \code{'GENENAME'}...)
#' @param display.identifier The name of the column of \code{annotation.db} used to label genes/features in \code{labelled.genes}
#' @export
volcano <- function(res.df,
                    alpha=0.1,
                    labelled.genes=NULL,
                    labelled.genes.keytype=NULL,
                    annotation.db=NULL,
                    identifier=NULL,
                    gene.info.columns=NULL,
                    display.identifier=NULL
                    ){
  args.list <- as.list(environment(), all=TRUE)

  # Get the annotation database and the identifier
  annotation.info <- annotation.db.helper(args.list,obj=res.df)

  if(is.null(labelled.genes.keytype)) labelled.genes.keytype <- annotation.info$identifier
  if(is.null(display.identifier)) display.identifier <- labelled.genes.keytype

  plot.colors <- with(res.df,ifelse(!is.null(padj) & padj < alpha, 'red','black'))
  res.df$log10p <- with(res.df,-1.0*log10(pvalue))

  aes_string_args_list <- list(
    "x"="log2FoldChange",
    "y"="log10p",
    "gene"=display.identifier,
    "log2FoldChange"="log2FoldChange",
    "baseMean"="baseMean",
    "padj"="padj"
  )
  info.columns <- as.list(annotation.info$gene.info.columns)
  names(info.columns) <- info.columns
  aes_string_args_list <- c(info.columns,aes_string_args_list)
  aes_string_args_list <- aes_string_args_list[!duplicated(aes_string_args_list)]

  g <- ggplot2::ggplot(res.df,
              do.call(ggplot2::aes_string,aes_string_args_list)
              ) +
                ggplot2::geom_point(color=plot.colors) +
                ggplot2::xlab(expression(paste(log[2],FC))) + ggplot2::ylab(expression(paste(-log[10],p))) +
                ggplot2::theme(axis.text=element_text(size=14), axis.title=ggplot2::element_text(size=16,face="bold"))
  if (!is.null(labelled.genes)){
    res.labelled <- res.df[labelled.genes,]
    g <- g + ggrepel::geom_text_repel(data=res.labelled,ggplot2::aes_string(label=display.identifier),color='red')
  }
  g
}
