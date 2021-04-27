#' Derive (guess) sample data from sample file names
#'
#' @param filenames Vector of sample file names.
#' @param splitchars Vector of delimiter characters for splitting filenames
#' @export
guess.sample.data <- function(filenames=NULL,splitchars=c('.','_','-')){
  if(is.null(filenames)){
    counts <- read.table('data/gene_counts.tsv',sep='\t',header=TRUE,row.names=1)
    filenames <- colnames(counts)
  }
  split.pattern <- paste(c('[',splitchars,']'),collapse='')
  df <- t(as.data.frame(strsplit(filenames,split=split.pattern)))
  df <- as.data.frame(df)
  rownames(df) <- filenames
  cols <- 1:ncol(df)
  #Is a particular column informative (i.e. does it vary over the samples)?
  is.informative <- function(i) length(unique(df[,i])) > 1
  df <- df[,sapply(cols,is.informative),drop=FALSE]
  rownames(df) <- apply(df, 1, function(l) paste(l,collapse='.'))
  df
}
