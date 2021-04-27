#' Do a complete comparison of two groups - differential expression, GO enrichment, and KEGG enrichment
#'
#' @param dds The DESeq2 object to use for the analysis
#' @param attr The attribute to be compared (i.e. 'genotype' or 'treatment')
#' @param level1 The value of the comparison attribute in group 1 (i.e. 'Knockout' or 'Infected')
#' @param level2 The value of the comparison attribute in group 1 (i.e. 'Wildtype' or 'Healthy')
#' @param output_directory The path to the output directory for the results
#' @param overwrite=FALSE Whether to overwrite previous results if output_directory already exists
#' @param annotation.db The annotation database (i.e. \code{org.Mm.eg.db}) used to look up information about features (genes or transcripts)
#' @param identifier The name of the column of \code{annotation.db} that corresponds to the row names of \code{dds} (feature names)
#' @param gene.info.columns Vector of column names from \code{annotation.db} to include in the output (for example, \code{'SYMBOL'}, \code{'GENENAME'}...)
#' @export
compare <- function(dds,
                        attr,level1,level2,
                        output_directory,
                        overwrite=FALSE,
                        annotation.db=NULL,
                        identifier=NULL,
                        gene.info.columns=NULL
                      ){
  #Get annotation information
  args.list <- as.list(environment(), all=TRUE)
  annotation.info <- annotation.db.helper(args.list,obj=dds)

  `%>%` <- magrittr::`%>%`

  alpha <- 0.1

  #Make sure that level1 and level2 are valid levels of the attribute
  for (l in c(level1,level2)){
    if (!(l %in% levels(colData(dds)[,attr]))){
      message <- paste0(sQuote(l)," is not a valid level of the sample data column ", sQuote(attr), ".")
      stop(message)
    }
  }

  if (file.exists(output_directory)){
    if(!overwrite){
      warning.message <- paste(c('Did not do the comparison of ', attr, ': ', level1, ' vs ', level2, ' because output directory already exists.'),sep='')
      stop(warning.message)
    } else {
      warning.message <- paste(c('Overwriting contents of ', output_directory, ' with the results of the comparison ',attr, ': ', level1, ' vs ', level2))
      warning(warning.message)
    }
  } else{
    dir.create(output_directory)
  }

  #Differential expression analysis
  de.path <- file.path(output_directory,'DifferentialExpression')
  if(!file.exists(de.path)) dir.create(de.path)

  de.up.file <- paste0(c('up_',level1,'_vs_',level2,'.tsv'),collapse='')
  de.down.file <- paste0(c('up_',level2,'_vs_',level1,'.tsv'),collapse='')

  reportDirectory <- file.path(de.path,paste0(c(level1,'_vs_',level2,'_HTML'),collapse=''))
  de.html.report(dds=dds,
                 output.dir = reportDirectory,
                 attr=attr,level1=level1,level2=level2,
                 annotation.db=annotation.info$annotation.db,
                 identifier=annotation.info$identifier,
                 gene.info.columns=annotation.info$gene.info.columns
  )

  res <- ggp.results(dds,
                     contrast=c(attr,level1,level2),
                     annotation.db=annotation.info$annotation.db,
                     identifier=annotation.info$identifier,
                     gene.info.columns=annotation.info$gene.info.columns
                    )
  res %>% dplyr::filter(!is.na(padj)) -> background
  background.file <- file.path(de.path,'background.tsv')
  write.table(background,background.file,
              sep='\t',
              quote=FALSE,
              row.names=FALSE)
  background %>% dplyr::filter(!is.na(padj) & padj < 0.1 & log2FoldChange > 0) %>% dplyr::arrange(-log2FoldChange) -> de.up
  background %>% dplyr::filter(!is.na(padj) & padj < 0.1 & log2FoldChange < 0) %>% dplyr::arrange(log2FoldChange) -> de.down
  write.table(de.up,file=file.path(de.path,de.up.file),sep = '\t',row.names = FALSE,quote=FALSE)
  write.table(de.down,file=file.path(de.path,de.down.file),sep = '\t',row.names = FALSE,quote=FALSE)

  #Volcano plot
  g.volcano <- volcano(res,
                       annotation.db=annotation.info$annotation.db,
                       identifier=annotation.info$identifier,
                       gene.info.columns=annotation.info$gene.info.columns,
                       display.identifier=annotation.info$display.identifier
                      )
  volcano.file.name <- file.path(output_directory,paste0(c('volcano_',level1,'_vs_',level2,'.pdf'),collapse=''))
  pdf(volcano.file.name)
  suppressWarnings(print(g.volcano))
  invisible(dev.off())

  #Interactive volcano plot
  #For some reason htmlwidgets::saveWidget only works with absolute paths
  output.directory.normalized <- normalizePath(output_directory)
  volcano.html.file.name <- file.path(output.directory.normalized,paste0(c('volcano_',level1,'_vs_',level2,'.html'),collapse=''))
  #Plotly can't handle subscripts, so redo the axis labels
  g.volcano.no.subscripts <- g.volcano + xlab('log_2 FC') + ylab("-log_10 p")
  volcano.widget <- plotly::ggplotly(g.volcano.no.subscripts)
  htmlwidgets::saveWidget(widget=volcano.widget,file=volcano.html.file.name,selfcontained=TRUE)

  #Want to label the top 10 DE genes (by logFC) in both direction
  #HOWEVER if there are fewer than 10 DE genes in either direction, be careful because
  #R is perfectly happy to draw 10 NAs from an empty data frame. See notes of 6/23/17
  #The head function will deal with this correctly
  top.genes.up <- head(de.up,10)$ENTREZID
  top.genes.down <- head(de.down,10)$ENTREZID
  top.20 <- c(top.genes.up,top.genes.down)
  scatter.plot.filename <- file.path(output_directory,paste0(c('scatterplot_',level1,'_vs_',level2,'.pdf'),collapse=''))
  g.scatter <- de.scatter.plot(dds,
                               attr,level1,level2,
                               labelled.genes = top.20,
                               annotation.db=annotation.info$annotation.db,
                               identifier=annotation.info$identifier,
                               gene.info.columns=annotation.info$gene.info.columns,
                               display.identifier=annotation.info$display.identifier
                              )
  pdf(scatter.plot.filename)
  print(g.scatter)
  invisible(dev.off())

  #Interactive DE scatter plot
  g.scatter.interactive <- de.scatter.plot(dds,
                               attr,level1,level2,
                               annotation.db=annotation.info$annotation.db,
                               identifier=annotation.info$identifier,
                               gene.info.columns=annotation.info$gene.info.columns,
                               display.identifier=annotation.info$display.identifier,
                               tick.label.superscripts = FALSE
  )
  scatter.html.file.name <- file.path(output.directory.normalized,paste0(c('scatterplot_',level1,'_vs_',level2,'.html'),collapse=''))
  scatter.widget <- plotly::ggplotly(g.scatter.interactive)
  htmlwidgets::saveWidget(widget=scatter.widget,file=scatter.html.file.name,selfcontained = TRUE)
}
