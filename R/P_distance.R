#' Percent difference Distance
#'
#' The P-distance is simply the percentage of the aligned base pairs where the
#' two bases are different.  See May 2004 for a discussion on the many ways to
#' compute percent identity (and therefore difference.)
#'
#' @param train_parameters A list containing the probabilities of the alignment,
#' produced by the [`readTrainFile()`] function.
#'
#' @returns A numeric value representing the evolutionary distance between two
#' genomes.  The larger the value, the more divergent the genomes.
#'
#' @family Similarity indexes
#'
#' @author Charles Plessy
#'
#' @references May, A. C. W. (2004). Percent Sequence Identity: The Need to Be Explicit.
#' *Structure*, 12(5), 737–738. \doi{doi.org/10.1016/j.str.2004.04.001}
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' P_distance(parameters)
#'
#' @export

P_distance <- function(train_parameters) {
  P <- train_parameters$probability_matrix
  1 - sum(diag(P))
}
