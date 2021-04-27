#' Import mapping rates from featureCounts report file
#'
#' @param report.file Path to the featureCounts report file (default: 'reports/featureCounts.report.txt')
#' @export
#'
import.mapping.rates <- function(report.file='reports/featureCounts.report.txt'){
  command.filenames <- paste0(
    'cat ',
    report.file,
    ' | grep "Process BAM file"',
    "| cut -d $' ' -f5 "
  )
  command.mapping.rate <- paste0(
    'cat ',
    report.file,
    ' | grep "Successfully assigned"',
    "| cut -d $' ' -f10"
  )
  filenames <- system(command.filenames,intern = TRUE)
  pieces <- strsplit(filenames,'/')
  filenames <- sapply(pieces,function(v) tail(v,n=1))
  filenames <- sapply(strsplit(filenames,".bam"),'[[',1)
  mapping.rate.strings <- system(command.mapping.rate,intern = TRUE)
  mapping.rates <- as.numeric(gsub("%\\)","",gsub("\\(","",mapping.rate.strings)))/100.0
  data.frame(
    filename=filenames,
    mapping.rate = mapping.rates
  )
}
