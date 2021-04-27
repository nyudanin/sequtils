#' Makes a GSEA plot from the output of the GSEA program. This is adapted from Thomas Kuilman's code, available at https://github.com/PeeperLab/Rtoolbox/blob/master/R/ReplotGSEA.R.
#'
#' @param path Path to the GSEA output directory
#' @param gene.set Name of the gene set. Note that this function uses a call to grep to find a matching gene set name.
#' @param main.title Main title of the plot (default: empty string).
#' @param show.p.value Whether to show the nominal p value (default: TRUE)
#' @param p.value.cutoff Lower cutoff below which p values should be displayed as inequalities only.
#' @param show.FDR Whether to show the False Discovery Rate (default: TRUE)
#' @param FDR.cutoff Lower cutoff below which FDRs should be displayed as inequalities only
#' @param balanced.color.scheme Whether to use a color scheme where 0 is white by definition (Default: FALSE)
#' @param left.label Label to be placed on the left side of the color gradient (default: 'Positive').
#' @param right.label Label to be placed on the right side of the color gradient (default: 'Negative').
#' @param show.metric.plot Boolean value: whether to show the lower plot, which shows the distribution of ranking metrics (default: TRUE).
#' @param ranking.metric.name Description of the ranking metric, for example 'Test statistic'
#' @export
replot.GSEA <- function(path,
                        gene.set,main.title='',
                        show.p.value=TRUE,
                        p.value.cutoff,
                        show.FDR=TRUE,
                        FDR.cutoff,
                        balanced.color.scheme=FALSE,
                        left.label='Positive',
                        right.label='Negative',
                        show.metric.plot=TRUE,
                        ranking.metric.name
                      ) {

  class.name=''

  if(missing(path)) {
    stop("Path argument is required")
  }
  if (!file.exists(path)) {
    stop("The path folder could not be found. Please change the path")
  }
  if(missing(gene.set)) {
    stop("Gene set argument is required")
  }

  ## Load .rnk data
  path.rnk <- list.files(path = file.path(path, "edb"),
                         pattern = ".rnk$", full.names = TRUE)
  gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")

  ## Load .edb data
  path.edb <- list.files(path = file.path(path, "edb"),
                         pattern = ".edb$", full.names = TRUE)
  gsea.edb <- read.delim(file = path.edb,
                         header = FALSE, stringsAsFactors = FALSE)
  gsea.edb <- unlist(gsea.edb)
  gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
  gsea.metric <- unlist(strsplit(gsea.metric, " "))
  gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
  gsea.metric <- gsub("METRIC=", "", gsea.metric)
  gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]

  # Select the right gene set
  if (length(gsea.edb) == 0) {
    stop(paste("The gene set name was not found, please provide",
               "a correct name"))
  }
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
    warning(paste("More than 1 gene set matched the gene.set",
                  "argument; the first match is plotted"))
  }
  gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]

  # Get template name
  gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
  gsea.edb <- unlist(strsplit(gsea.edb, " "))
  gsea.template <- gsea.edb[1]

  # Get gene set name
  gsea.gene.set <- gsea.edb[2]
  gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)

  # Get enrichment score
  gsea.enrichment.score <- gsea.edb[3]
  gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)

  # Get gene set name
  gsea.normalized.enrichment.score <- gsea.edb[4]
  gsea.normalized.enrichment.score <- gsub("NES=", "",
                                           gsea.normalized.enrichment.score)

  # Get nominal p-value
  gsea.p.value <- gsea.edb[5]
  gsea.p.value <- gsub("NP=", "", gsea.p.value)
  if(as.numeric(gsea.p.value) < p.value.cutoff){
    gsea.p.value.string <- paste0('< ',p.value.cutoff)
  } else{
    gsea.p.value.string <- paste0('= ',gsea.p.value)
  }

  p.value.display.string <- ifelse(show.p.value,paste0("Nominal p-value ", gsea.p.value.string,' \n'),'')

  # Get FDR
  gsea.fdr <- gsea.edb[6]
  gsea.fdr <- gsub("FDR=", "", gsea.fdr)
  gsea.fdr <- as.numeric(gsea.fdr)

  if(gsea.fdr < FDR.cutoff){
    fdr.string <- paste0('< ',FDR.cutoff)
  } else{
    fdr.string <- paste0('= ',gsea.fdr)
  }

  fdr.display.string <- ifelse(show.FDR,paste0("FDR ",fdr.string, ' \n'),'')

  # Get hit indices
  gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
  gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
  gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
  gsea.hit.indices <- as.integer(gsea.hit.indices)

  # Get ES profile
  gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
  gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
  gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
  gsea.es.profile <- as.numeric(gsea.es.profile)


  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)

  # Create a new device of appropriate size
  #dev.new(width = 3, height = 3)

  # Create a division of the device
  if (show.metric.plot){
    gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 1))
  } else{
    gsea.layout <- layout(matrix(c(1, 2, 3)), heights = c(1.7, 0.5, 0.2))
  }

  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
       c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",cex.lab=1.5,
       main = list(main.title, font = 1, cex = 1.5),
       panel.first = {
         abline(h = seq(round(min(gsea.es.profile), digits = 1),
                        max(gsea.es.profile), 0.1),
                col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
  plot.coordinates <- par("usr")
  if(gsea.enrichment.score < 0) {
    text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
         paste(p.value.display.string, fdr.display.string, "ES:",
               prettyNum(as.numeric(gsea.enrichment.score),digits=3), "\nNormalized ES:",
               prettyNum(as.numeric(gsea.normalized.enrichment.score),digits=3)), adj = c(0, 0),cex=1.2)
  } else {
    text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
         paste(p.value.display.string,fdr.display.string, "ES:",
               prettyNum(as.numeric(gsea.enrichment.score),digits=3), "\nNormalized ES:",
               prettyNum(as.numeric(gsea.normalized.enrichment.score),digits=3), "\n"), adj = c(1, 1),cex=1.2)
  }

  par(mar = c(0, 5, 0, 2))
  plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
       ylab = "", xlim = c(1, length(gsea.rnk$metric)))
  abline(v = gsea.hit.indices, lwd = 0.75)

  par(mar = c(0, 5, 0, 2))

  rank.colors <- gsea.rnk$metric
  if(balanced.color.scheme){
    rank.colors <- ifelse(rank.colors > 0, rank.colors/max(rank.colors), rank.colors/abs(min(rank.colors)))
  }

  rank.colors <- rank.colors - min(rank.colors)
  rank.colors <- rank.colors / max(rank.colors)
  rank.colors <- ceiling(rank.colors * 255 + 1)
  rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n",
          xlim = c(1, length(gsea.rnk$metric))
  )
  box()
  text(length(gsea.rnk$metric) / 2, 0.7,
       labels = ifelse(!missing(class.name), class.name, gsea.template))
  text(length(gsea.rnk$metric) * 0.01, 0.7, left.label, adj = c(0, NA),cex=1.2)
  text(length(gsea.rnk$metric) * 0.99, 0.7, right.label, adj = c(1, NA),cex=1.2)

  if(show.metric.plot){

    par(mar = c(5, 5, 0, 2))
    rank.metric <- rle(round(gsea.rnk$metric, digits = 2))

    metric.min <- min(gsea.rnk$metric)
    metric.max <- max(gsea.rnk$metric)

    plot(gsea.rnk$metric, type = "n", xaxs = "i",
         xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
         ylim = c(metric.min, metric.max), yaxs = "i",
         ylab = '',#if(gsea.metric == "None") {ranking.metric.name} else {gsea.metric},
         cex.lab=1.5
         #panel.first = abline(h = seq(-0.5, 0.5, 0.5), col = "red", lty = 2)
    )

    barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
            xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
            ylim = c(metric.min, metric.max), yaxs = "i", width = rank.metric$lengths, border = NA,
            ylab = ifelse(gsea.metric == "None", ranking.metric.name, gsea.metric), space = 0, add = TRUE,
            cex.lab=1.5
    )
    box()
  }

  # Reset to default
  #par(def.par)
}
