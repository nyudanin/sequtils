#' Scatter plot of normalized counts
#'
#' Makes a scatter plot comparing the mean normalized counts in two groups of samples. Alternatively if the parameter \code{returnData=TRUE} is
#' supplied, returns the dataframe used to make the plot. Genes with significantly higher expression in the group where \code{attr==level1} are
#' highlighted in red, while genes with significantly higher expression in the group where \code{attr=level2} are highlighted in red. Optionally,
#' a list of genes \code{labelled.genes} may be given, which which be labelled in the plot. The labelled genes will be looked up in the annotation
#' database \code{annotation.db} using the column named \code{labelled.genes.keytype}. As in many \code{ggp.rnaseq} functions, \code{annotation.db}
#' can be supplied via the \code{'annotation.db'} attribute of the object \code{dds}, and will default to \code{org.Mm.eg.db}. Once the genes
#' \code{labelled.genes} have been selected, they are labelled in the plot using the column \code{display.identifier} of \code{annotation.db}.
#' If \code{display.identifier} is not supplied, \code{labelled.genes.keytype} is used.
#'
#' @param dds DESeq object
#' @param attr The attribute to be compared (i.e. 'genotype' or 'treatment')
#' @param level1 The value of the comparison attribute in group 1, which will be on the y axis (i.e. 'Knockout' or 'Infected')
#' @param level2 The value of the comparison attribute in group 2, which will be on the x axis (i.e. 'Wildtype' or 'Healthy')
#' @param alpha Significance cutoff for highlighting genes
#' @param labelled.genes Vector of identifiers of genes to be labelled
#' @param labelled.genes.keytype The column name (e.g. \code{'SYMBOL'}) from \code{annotation.db} to be used to look up \code{labelled.genes}
#' @param annotation.db The annotation database (i.e. \code{org.Mm.eg.db}) used to look up information about features (genes or transcripts)
#' @param gene.info.columns Vector of column names from \code{annotation.db} to include in the output (for example, \code{'SYMBOL'}, \code{'GENENAME'}...)
#' @param identifier The name of the column of \code{annotation.db} that corresponds to the row names of \code{dds} (feature names)
#' @param display.identifier The name of the column of \code{annotation.db} used to label genes/features in \code{labelled.genes}
#' @param pseudocount Pseudocount to add to all gene/feature counts before taking the logarithm (default: 0.5)
#' @param label.color Color with which to label genes (default: \code{'orange'})
#' @param returnData Boolean (FALSE by default). If TRUE, return the dataframe used to make the plot, rather than the plot itself
#' @param tick.label.superscripts Boolean (TRUE by default) for whether to render tick labels with superscripts (these don't render in interactive plots!)
#' @export
de.scatter.plot <- function(dds,
                            attr,
                            level1,
                            level2,
                            alpha=0.1,
                            labelled.genes=NULL,
                            labelled.genes.keytype=NULL,
                            annotation.db=NULL,
                            identifier=NULL,
                            gene.info.columns=NULL,
                            display.identifier=NULL,
                            pseudocount=0.5,
                            label.color='orange',
                            returnData = FALSE,
                            tick.label.superscripts=TRUE
                          ){
  args.list <- as.list(environment(), all=TRUE)

  # First, get the annotation database and the identifier
  annotation.info <- annotation.db.helper(args.list,obj=dds)

  if(is.null(labelled.genes.keytype)) labelled.genes.keytype <- annotation.info$identifier
  if(is.null(display.identifier)) display.identifier <- labelled.genes.keytype

  #gene.info.columns probably may not be set in either the args.list or the attributes of dds
  #In case it's not defined, just make sure to include the identifier and the display identifier
  if(is.null(annotation.info$gene.info.columns)){
    annotation.info$gene.info.columns <- unique(c(annotation.info$identifier,display.identifier))
  }

  res.df <- ggp.results(dds,
                        contrast=c(attr,level1,level2),
                        annotation.db = annotation.info$annotation.db,
                        identifier=annotation.info$identifier,
                        gene.info.columns = annotation.info$gene.info.columns)

  counts.df <- counts(dds,normalize=TRUE)
  groupMeans <- cbind(
    rowMeans(counts.df[,colData(dds)[,attr] == level1]) + pseudocount, #Add a psuedocount before we take logs
    rowMeans(counts.df[,colData(dds)[,attr] == level2] + pseudocount)
  )
  colnames(groupMeans) <- paste("norm.counts",c(level1,level2),sep = '.')
  df <- data.frame(res.df,groupMeans)

  if(returnData) return(df)

  plot.colors <- rep('black',nrow(df))
  plot.colors[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange > 0] <- 'red'
  plot.colors[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange < 0] <- 'blue'
  plot.alphas <- ifelse(!is.na(df$padj) & df$padj < alpha,1.0,0.1)

  aes_string_args_list <- list(
    "x"=paste0('norm.counts.',level2),
    "y"=paste0('norm.counts.',level1),
    "label"=display.identifier,
    "log2FoldChange"="log2FoldChange",
    "baseMean"="baseMean",
    "padj"="padj"
  )
  info.columns <- as.list(annotation.info$gene.info.columns)
  names(info.columns) <- info.columns
  aes_string_args_list <- c(info.columns,aes_string_args_list)
  aes_string_args_list <- aes_string_args_list[!duplicated(aes_string_args_list)]
  g <- ggplot2::ggplot(
    df,
    do.call(ggplot2::aes_string,aes_string_args_list)
  ) + ggplot2::geom_point(color=plot.colors,alpha=plot.alphas) +
                            ggplot2::xlab(paste0("Mean normalized counts: ",level2)) +
                            ggplot2::ylab(paste0("Mean normalized counts: ",level1)) +
                            ggplot2::theme(axis.text=element_text(size=12),axis.title=ggplot2::element_text(size=14)
                        )

  if(tick.label.superscripts){
    g <- g + ggplot2::scale_x_continuous(trans='log10',labels=scales::trans_format('log10',scales::math_format(10^.x))) +
             ggplot2::scale_y_continuous(trans='log10',labels=scales::trans_format('log10',scales::math_format(10^.x)))
  } else{
    g <- g + ggplot2::scale_x_log10() + ggplot2::scale_y_log10()
  }

  if(length(labelled.genes)==0){
    labelled.genes <- NULL
  }
  if(!is.null(labelled.genes)){
    bad.genes <- setdiff(labelled.genes,df[,labelled.genes.keytype])
    if(length(bad.genes) > 0){
      warning(paste0("The following genes were not found by de.scatter.plot: ",paste(as.character(bad.genes),collapse=', ')))
    }
    labelled.genes.good <- intersect(df[,labelled.genes.keytype],labelled.genes)
    df.labelled <- df[df[,labelled.genes.keytype] %in% labelled.genes.good,]
    g <- g +
      ggrepel::geom_text_repel(
      data=df.labelled,
      ggplot2::aes_string(
        x=paste0('norm.counts.',level2),
        y=paste0('norm.counts.',level1),
        label=display.identifier
      ),
      force=20,
      color=label.color
    )
  }
  g
}
