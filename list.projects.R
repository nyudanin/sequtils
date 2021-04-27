#' List project data packages
#'
#' @export
list.projects <- function(){
  package.names <- as.data.frame(installed.packages())$Package
  as.character(package.names[grep("project",package.names)])
}
