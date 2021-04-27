#' Dummy plot for when a gene symbol or Entrez ID is not found in the database
#'
#' @param identifier Identifier (symbol or Entrez ID) of the gene
gene.not.found.plot <- function(identifier){
  df <- data.frame(x=c(1,1),y=c(1.5,1),text=c("Gene not found:",identifier))
  ggplot(df,aes(x,y,label=text)) + geom_text(size=20) + scale_x_continuous(limits=c(0,2)) + scale_y_continuous(limits=c(0,2))
}

#' Plot counts for a given gene
#'
#' @param entrez Entrez ID of gene
#' @param symbol Gene symbol
#' @param normalized Boolean: whether or not to give the normalized counts rather than the raw counts
#' @param labelled Boolean: whether to label the samples
#' @param sample.data Data frame of sample data
#' @param x Column name from sample.data to be used for the x axis
#' @param color Column name from sample.data to be used for the color
#' @param color.scheme Named list giving the color mapping
#' @export
gene.plot.all <- function(entrez=NULL,
                          symbol=NULL,
                          normalized=TRUE,
                          labelled=FALSE,
                          sample.data=coldata,
                          x='treatment',
                          color='treatment',
                          color.scheme=NULL
                          ){
  suppressPackageStartupMessages(OK <- require(ggplot2))
  if (!OK) stop("Error: ggplot2 package not found")
  suppressPackageStartupMessages(OK <- require(org.Hs.eg.db))
  if (!OK) stop("Error: org.Hs.eg.db package not found")
  suppressPackageStartupMessages(OK <- require(DESeq2))
  if (!OK) stop("Error: DESeq2 package not found")

  if (!is.null(entrez)){
    if (!(symbol %in% AnnotationDbi::keys(org.Mm.eg.db,keytype = 'ENTREZID'))){
      return(gene.not.found.plot(symbol))
    }
    suppressMessages(
      symbol <- AnnotationDbi::select(org.Mm.eg.db,entrez,keytype="ENTREZID",columns='SYMBOL')$SYMBOL
    )
  }
  else{
    if (!is.null(symbol)){
      if (!(symbol %in% AnnotationDbi::keys(org.Mm.eg.db,keytype = 'SYMBOL'))){
        return(gene.not.found.plot(symbol))
      }
      suppressMessages(
        entrez <- AnnotationDbi::select(org.Mm.eg.db,symbol,keytype="SYMBOL",columns='ENTREZID')$ENTREZID
      )
    }
    else{
      stop("Need either an Entrez ID or a gene symbol as input")
    }
  }

  if(all(as.numeric(countdata[entrez,]) == 0)){
    df <- data.frame(
      sample.data,
      counts=0,
      label=rownames(sample.data)
    )
    g <- ggplot(df,aes_string(x=x,y='counts',label='label',color=color)) + scale_y_continuous(limits=c(0,100))
    if (labelled){
      g <- g + geom_text()
    }
    g <- g + ggtitle(symbol) + labs("Raw counts")
    return(g)
  }

  df <- data.frame(
    sample.data,
    counts=DESeq2::counts(dds,normalized=normalized)[entrez,rownames(sample.data)],
    label=rownames(sample.data)
  )
  # gene.info <- info.from.entrez(entrez)
   normString <- ifelse(normalized,"Normalized","Raw")
   m <- max(df$counts)
   g <- ggplot(df,aes_string(x=x,y='counts',label='label',color=color))
   if (labelled){
    g <- g + geom_text()
  }
  else{
    g <- g + geom_jitter(height=0,width=0.1,size=3) #+
      #scale_fill_manual(values=treatmentColors)
  }
  if (!is.null(color.scheme)){
    color.values <- color.scheme[sample.data[,x]]
    #g <- g + scale_fill_manual(values=color.values)
    g <- g + scale_color_manual(values=color.values)
  }
  g <- g + ggtitle(symbol) + labs(y=paste(normString,'Counts')) + scale_y_continuous(limits=c(0,1.1*m))
  g
}
