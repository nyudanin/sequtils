#' Make a plot of the CDF (cumulative distribution function) of gene expression levels
#'
#' @param countdata Data frame of expression levels (counts). The rows should be genes/transcripts and the columns should be samples
#' @export
expression.cdf <- function(countdata){
  suppressPackageStartupMessages(OK <- require(ggplot2))
  if (!OK) stop("Error: ggplot2 package not found")

}
