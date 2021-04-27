#' Scatterplot matrix
#'
#' @param dds DESeq2 object (post-normalization)
#' @param coldata Sample data frame
#' @param group Name of column of interest in sample data frame
#' @param output.dest Path to directory for output plots
#' @param highlight.genes List of Entrez IDs of genes to be highlighted
#' @export
scatterplot.matrix <- function(dds,coldata,group,output.dest,highlight.genes=c()){
  suppressPackageStartupMessages(OK <- require(DESeq2))
  if (!OK) stop("Error: package 'DESeq2' not found")
  suppressPackageStartupMessages(OK <- require(ggplot2))
  if (!OK) stop("Error: package 'ggplot2' not found")
  suppressPackageStartupMessages(OK <- require(reshape2))
  if (!OK) stop("Error: package 'reshape2' not found")

normalized.counts <- as.data.frame(counts(dds,normalize=TRUE))
  highlight.genes <- as.character(highlight.genes)

  groups <- unique(coldata[,group])
  for (i in seq(length(groups))){
    for (j in seq(i,length(groups))){
      group1 <- groups[i]
      group2 <- groups[j]
      subjects.1 <- rownames(coldata)[coldata[,group] == group1]
      subjects.2 <- rownames(coldata)[coldata[,group] == group2]
      temp <- normalized.counts[,subjects.1]
      temp$entrez <- rownames(temp)
      counts.long <- melt(temp,variable.name='subject.1',value.name='count.1',
                          varnames=subjects.1,id.vars='entrez'
      )
      counts.long <- cbind(counts.long,normalized.counts[counts.long$entrez,subjects.2])
      counts.long.2 <- melt(counts.long,variable.name='subject.2',value.name='count.2',
                            varnames=subjects.2,id.vars=c('entrez','subject.1','count.1')
      )
      plot.name <- file.path(output.dest,
                             paste0(
                               paste(group1,group2,sep = '_'),
                               '.png'
                               )
                            )

      plot.colors <- ifelse(counts.long.2$entrez %in% highlight.genes,'red','black')

      png(plot.name,width=7,height=7,units='in',res=150)
      g <- ggplot(
        counts.long.2,
        aes(log(count.1),log(count.2))
      ) + geom_point(color=plot.colors) +
        geom_abline(slope=1,color='red',linetype='dashed') +
        facet_grid(subject.2 ~ subject.1)
      print(g)
      invisible(dev.off())
    }
  }
}
