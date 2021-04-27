#' Create a bash script to run GSEA (Preranked mode)
#'
#' @param gsea.jar.path Path to the GSEA executable .jar file. Default: '~/programs/gsea/gsea2-2.2.2.jar'
#' @param gmt.file Path to GMT file (gene set database)
#' @param ranked.gene.file Path to RNK file (ranked gene list)
#' @param n.perm Number of permutations (default: 1000)
#' @param scoring.scheme Scoring scheme for the enrichment score. Choices are 'classic','weighted','weighted_p2', or 'weighted_p1.5'. Default: 'classic'.
#' @param output.prefix Filename prefix for output
#' @param plot.top.n Number of the top gene sets for which to make plots (default: 20).
#' @param max.size Maximum size of gene sets to include in the analysis (default: 5000).
#' @param min.size Minimum size of gene sets to include in the analysis (default: 5).
#' @param filename File name of the output bash script
#' @export
write.gsea.script.ranked <- function(
  gsea.jar.path='~/programs/gsea/gsea2-2.2.2.jar',
  gmt.file,
  ranked.gene.file,
  n.perm=1000,
  scoring.scheme='classic',
  output.prefix,
  plot.top.n=20,
  max.size=5000,
  min.size=5,
  filename
){
  command.part1 <- sprintf("java -Xmx4096m -cp %s xtools.gsea.GseaPreranked -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm %s \\",
                     gsea.jar.path,
                     gmt.file,
                     n.perm
                     )
  command.part2 <- sprintf("        -rnk %s -scoring_scheme %s -rpt_label %s -plot_top_x %s -rnd_seed timestamp -set_max %s -set_min %s -zip_report false -out . -gui false",
                          ranked.gene.file,
                          scoring.scheme,
                          output.prefix,
                          plot.top.n,
                          max.size,
                          min.size
                          )
  command <- paste(command.part1,command.part2,collapse=' ')
  #cat(command,'\n')
  fileConn<-file(filename)
  writeLines(c('#!/bin/bash',command.part1,command.part2), fileConn)
  close(fileConn)

  system(paste0('chmod +x ',filename))
}
