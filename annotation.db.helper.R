#' Internal function to determine the annotation database (\code{annotation.db}), gene identifier column (\code{identifier}),
#' and gene info columns (\code{gene.info.columns}).
#'
#' This function will first look in \code{args.list} for this information.
#' If it doesn't find \code{annotation.db} in \code{args.list}, it will seek all of this information in the attributes of
#' \code{obj} (usually either a \code{DESeqDataSet} or \code{DESeqResults}). Finally, if there is no attribute \code{annotation.db}
#' for the object \code{obj}, it will fall back on default values, namely \code{annotation.db = org.Mm.eg.db}, \code{identifier = 'ENTREZID'},
#' and \code{gene.info.columns=c('ENTREZID','SYMBOL','GENENAME')}
#' @param args.list A list of arguments - highest priority source of values for \code{annotation.db}, \code{identifier}, and \code{gene.info.columns}
#' @param obj An object whose attributes are used to find values of \code{annotation.db}, etc., if \code{annotation.db} is not in \code{args.list}
annotation.db.helper <- function(args.list, obj){
  annotation.info <- list()
  if('annotation.db' %in% names(args.list) & !is.null(args.list[['annotation.db']])){
    annotation.info[['annotation.db']] <- args.list[['annotation.db']]

    #If we get annotation.db from the arguments list, we must get the other items as well
    if('identifier' %in% names(args.list) & !is.null(args.list[['identifier']])){
      annotation.info[['identifier']] <- args.list[['identifier']]
    }
    if('gene.info.columns' %in% names(args.list) & !is.null(args.list[['gene.info.columns']])){
      annotation.info[['gene.info.columns']] <- args.list[['gene.info.columns']]
    }
    if('display.identifier' %in% names(args.list) & !is.null(args.list[['display.identifier']])){
      annotation.info[['display.identifier']] <- args.list[['display.identifier']]
    }
    return(annotation.info)
  }

  #If we get here, annotation.db was not in the arguments list
  if('annotation.db' %in% names(attributes(obj))){
    annotation.info[['annotation.db']] <- attr(obj,'annotation.db')

    #If we get annotation.db from the attributes of the object, we get everything else there also
    if('identifier' %in% names(attributes(obj))){
      annotation.info[['identifier']] <- attr(obj,'identifier')
    }
    if('gene.info.columns' %in% names(attributes(obj))){
      annotation.info[['gene.info.columns']] <- attr(obj,'gene.info.columns')
    }
    if('display.identifier' %in% names(attributes(obj))){
      annotation.info[['display.identifier']] <- attr(obj,'display.identifier')
    }
    return(annotation.info)
  }

  #If we get here, we're falling back on default values
  suppressPackageStartupMessages(OK <- require(org.Hs.eg.db))
  if (!OK) stop("Error: org.Hs.eg.db package not found")
  annotation.info[['annotation.db']] <- org.Hs.eg.db
  annotation.info[['identifier']] <- 'ENTREZID'
  annotation.info[['gene.info.columns']] <- c('ENTREZID','SYMBOL','GENENAME')
  annotation.info[['display.identifier']] <- 'SYMBOL'

  #Make sure everything is of the correct type
  if(!(class(annotation.info[['annotation.db']]) %in% c('OrgDb','EnsDb'))){
    stop("Annotation database must be an OrgDb or EnsDb object")
  }
  if( class(annotation.info[['identifier']]) != 'character' ) stop("Gene identifier must be a string")
  if( class(annotation.info[['gene.info.columns']]) != 'character' ) stop("Gene info columns must be a vector of strings")

  #Make sure that identifier is among the columns of annotation.db
  if(!(annotation.info[['identifier']] %in% columns(annotation.info[['annotation.db']]))) stop("Gene identifier must be among columns of annotation database")
  if(!(all(annotation.info[['gene.info.columns']] %in% columns(annotation.info[['annotation.db']])))){
    stop("Gene info columns must be among the columns of the annotation database")
  }

  annotation.info
}
