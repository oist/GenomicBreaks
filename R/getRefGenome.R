#' @name getRefGenome
#'
#' @title Attempt to load a BSgenome
#' @details Internal function that retreives a BSgenome object or throws an error if not available. Adapted from the CAGEr package.
#' @return A BSgenome object
#' @param reference.genome name of genome package (and also genome data object)
#' @author Charles Plessy
#' @importFrom utils installed.packages
#' @export
#' @noRd


getRefGenome <- function(reference.genome) {

  if (is.null(reference.genome))

    stop("Can not run this function with a NULL genome; see ", sQuote('help("genomeName")'), ".")

  if(reference.genome %in% rownames(installed.packages()) == FALSE)

    stop("Requested genome is not installed! Please install required BSgenome package before running GenomicBreaks.")

  requireNamespace(reference.genome)

  getExportedValue(reference.genome, reference.genome)
}
