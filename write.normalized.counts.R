#' Write a table of normalized counts to a file
#'
#' @param dds Input DESeq2 object
#' @param filename Name of output file
#' @param sample.attributes Vector of names of sample attributes to be included (Default: none)
#' @param identifier The type of gene identifier ('ENTREZID','ENSEMBL'). Default is 'ENTREZID'.
#' @param overwrite Boolean: whether or not to overwrite a preexisting file. Default: FALSE
#' @export
write.normalized.counts <-
  function(dds,
           filename = 'data/normalized_counts.tsv',
           sample.attributes = NULL,
           identifier = 'ENTREZID',
           overwrite = FALSE) {
    suppressPackageStartupMessages(OK <- require(DESeq2))
    if (!OK)
      stop("Error: DESeq2 package not found")
    suppressPackageStartupMessages(OK <- require(org.Hs.eg.db))
    if (!OK)
      stop("Error: org.Hs.eg.db package not found")
    suppressPackageStartupMessages(OK <- require(dplyr))
    if (!OK)
      stop("Error: dplyr package not found")

    if (file.exists(filename)) {
      if (!overwrite) {
        warning.message <- paste(
          c(
            'Did not write to file ',
            filename,
            " because it exists already. Set overwrite=TRUE if you're sure you want to overwrite."
          ),
          sep = ''
        )
        stop(warning.message)
      } else {
        warning.message <- paste(c('Overwriting contents of ', filename))
        warning(warning.message)
      }
    }

    norm.counts.df <- as.data.frame(counts(dds, normalize = TRUE))

    if (identifier == 'ENTREZID') {
      entrez.ids <- rownames(norm.counts.df)
    }
    if (identifier == 'ENSEMBL') {
      entrez.ids <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = rownames(norm.counts.df),
        keytype = 'ENSEMBL',
        column = "ENTREZID",
        multiVals = 'first'
      )
      #This line would replace one of the samples' counts with ENSEMBL IDs!
      #norm.counts.df <- cbind(ENSEMBL=rownames(norm.counts.df),norm.counts.df)
    }

    gene.df <- suppressMessages(info.from.entrez(entrez.ids))

    if (!is.null(sample.attributes)) {
      sample.data.df <-
        as.data.frame(colData(dds))[, sample.attributes, drop = FALSE]
      sample.data.df <- mutate_all(sample.data.df, as.character)
      sample.order <- do.call(order, sample.data.df)
      sample.data.reordered <- sample.data.df[sample.order, ]
      norm.counts.df <- norm.counts.df[, sample.order]
      sample.data.t <- t(sample.data.reordered)
      sample.data.plus.counts <- rbind(sample.data.t,
                                       sapply(norm.counts.df, as.character))

      padding <- as.data.frame(matrix('', nrow = length(sample.attributes), ncol =
                                        3))

      #gene.df.plus.padding <- rbind(
      #  d1=c('','',''),
      #  d2=c('','',''),
      #  gene.df
      #)

      names(padding) <- colnames(gene.df)

      gene.df.plus.padding <- rbind(padding,
                                    gene.df)

      colnames(gene.df.plus.padding) <-
        c('ENTREZID', 'SYMBOL', 'GENENAME')
      gene.df.plus.padding <-
        mutate_all(gene.df.plus.padding, as.character)
      gene.df.plus.padding[1:length(sample.attributes), 'ENTREZID'] <-
        paste0(sample.attributes, ' -->')

      res <-
        cbind(gene.df.plus.padding,
              sample.data.plus.counts,
              row.names = NULL)

      norm.counts.df <- res
    }
    write.table(
      norm.counts.df,
      file = filename,
      sep = '\t',
      row.names = FALSE,
      quote = FALSE
    )
  }
