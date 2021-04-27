#' Make PCA plot or biplot
#'
#' @param x Matrix of input data: rows are genes and columns are samples
#' @param coldata Data frame of sample data
#' @param color.attr Name of column from coldata used to color samples
#' @param color.scheme Named list giving mapping of color.attr to colors
#' @param n.topvar.genes Number of genes (with highest variance) used to do PCA
#' @param biplot Boolean: whether or not to make a biplot
#' @param arrow.length Arbitrary length factor of biplot arrows
#' @param highlight.genes Vector of (character) Entrez gene IDs whose arrows will be highlighted
#' @param highlight.color Color with which to highlight genes from highlight.genes vector
#' @export
pca <- function(x,
                coldata,
                color.attr = NULL,
                color.scheme = NULL,
                n.topvar.genes = 500,
                biplot = FALSE,
                arrow.length = 0.9,
                highlight.genes = c(),
                highlight.color = 'red') {
  suppressPackageStartupMessages(OK <- require(genefilter))
  if (!OK)
    stop("Error: package 'genefilter' not found")
  suppressPackageStartupMessages(OK <- require(ggplot2))
  if (!OK)
    stop("Error: package 'ggplot2' not found")

  rv <- rowVars(assay(x))
  select <- order(rv, decreasing = TRUE)[1:n.topvar.genes]
  pca <- prcomp(t(assay(x)[select,]), center = TRUE)
  pc.data <- data.frame(obsnames = row.names(pca$x), pca$x, coldata)
  percent.variance <- round(100 * pca$sdev ** 2 / (sum(pca$sdev ** 2)))
  my.data <- data.frame(obsnames = row.names(pca$x), pca$x, coldata)
  datapc <- data.frame(varnames = rownames(pca$rotation), pca$rotation)
  mult <- min((max(my.data[, 'PC2']) - min(my.data[, 'PC2']) / (max(datapc[, 'PC2']) -
                                                                  min(datapc[, 'PC2']))),
              (max(my.data[, 'PC1']) - min(my.data[, 'PC1']) / (max(datapc[, 'PC1']) -
                                                                  min(datapc[, 'PC1']))))
  datapc <- transform(datapc,
                      v1 = arrow.length * mult * PC1,
                      v2 = arrow.length * mult * PC2)
  datapc$v <- with(datapc, sqrt(v1 * v1 + v2 * v2))
  g <- ggplot() +
    geom_point(data = pc.data,
               aes_string(x = 'PC1',
                          y = 'PC2',
                          color = color.attr))
  if (!is.null(color.scheme)) {
    g <- g + scale_color_manual(values = color.scheme)
  }
  if (biplot) {
    arrow.colors <-
      ifelse(rownames(datapc) %in% highlight.genes,
             highlight.color,
             "black")
    arrow.alphas <-
      ifelse(rownames(datapc) %in% highlight.genes, 1.0, 0.2)
    g <- g + geom_segment(
      data = datapc,
      aes(
        x = 0,
        y = 0,
        xend = v1,
        yend = v2
      ),
      arrow = arrow(length = unit(0.2, "cm")),
      color = arrow.colors,
      alpha = arrow.alphas
    )
  }
  g <- g + xlab(paste0("PC1: ", percent.variance[1], "% variance")) +
    ylab(paste0("PC2: ", percent.variance[2], "% variance"))
  print(g)
  invisible(list(
    gg = g,
    pc.data = pc.data,
    weights = pca$rotation
  ))
}
