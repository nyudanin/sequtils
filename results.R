#' Process a \code{DESeq2} results object into an annotated data frame
#'
#' This is a wrapper function for \code{DESeq2}'s results() function, meant to include additional information about genes/features.
#' The annotation database (either an \code{EnsDb} or an \code{OrgDb} object) must be specified. The \code{ggp.results} function will first
#' look for an `annotation.db` parameter. If that is not specified, it will look for an attribute of \code{dds} named \code{"annotation.db"}.
#' If that attribute is missing, it will fall back on `org.Mm.eg.db` as a default.
#' The parameter \code{identifier} should specify the name of the column in \code{annotation.db} corresponding to the feature/gene names -
#' that is, the rownames of \code{dds}. If \code{annotation.db} was specified using the attribute of \code{dds}, then \code{identifier} will
#' be the \code{"identifier"} attribute of \code{dds}.
#' The parameter \code{gene.info.columns} should be a vector of column names from `annotation.db`, specifying which information about
#' genes/features should be included in the output. If \code{annoation.db} was supplied using an attribute of \code{dds}, then
#' \code{gene.info.columns} will likewise be specified by the attribute \code{"gene.info.columns"} of \code{dds}.
#' The default fallback behavior is to use \code{annotation.db = org.Mm.eg.db}, \code{identifier='ENTREZID'}, and
#' \code{gene.info.columns=c('ENTREZID','SYMBOL','GENENAME')}
#' @param dds DESeq2 object, ready to be passed to DESeq2's results() function
#' @param contrast Object specifying the comparison, to be passed on to \code{DESeq2::results()}
#' @param annotation.db The annotation database (i.e. \code{org.Mm.eg.db}) used to look up information about features (genes or transcripts)
#' @param identifier The name of the column of \code{annotation.db} that corresponds to the row names of \code{dds} (feature names)
#' @param gene.info.columns Vector of column names from \code{annotation.db} to include in the output (for example, \code{'SYMBOL'}, \code{'GENENAME'}...)
#'
#' @export
ggp.results <- function(dds,contrast,annotation.db=NULL, identifier=NULL, gene.info.columns=NULL){
  args.list <- as.list(environment(), all=TRUE)

  if(!is(dds,"DESeqDataSet")) stop("Parameter dds must be a DESeq2 data set of class DESeqDataSet.")

  res <- DESeq2::results(dds,contrast)
  res.df <- as.data.frame(res)
  rownames(res.df) <- rownames(res)

  annotation.info <- annotation.db.helper(args.list,obj=dds)

  if(!(class(annotation.info$annotation.db) %in% c('OrgDb','EnsDb'))){
    stop("Parameter annotation.db must be an OrgDb or EnsDb object")
  }

  gene.info <- suppressMessages(AnnotationDbi::select(annotation.info$annotation.db,
                                                      keys=rownames(res.df),
                                                      keytype=annotation.info$identifier,
                                                      columns=annotation.info$gene.info.columns)
                                )

  indices <- match(rownames(res.df),gene.info[,annotation.info$identifier])

  #Return results data frame
  final.df <- data.frame(row.names=rownames(res.df),gene.info[indices,,drop=FALSE],res.df)

  #Pass on the annotation info as attributes of the resulting data frame
  attr(final.df,'annotation.df') <- annotation.info$annotation.db
  attr(final.df,'identifier') <- annotation.info$identifier
  attr(final.df,'gene.info.columns') <- annotation.info$gene.info.columns

  #Pass on the metadata from the DESeqResults object
  attr(final.df,'metadata') <- metadata(res)

  final.df
}
