#' Percent difference Distance
#'
#' The P-distance is simply the percentage of the aligned base pairs where the
#' two bases are different.  See May 2004 for a discussion on the many ways to
#' compute percent identity (and therefore difference.)
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @returns A numeric value representing the evolutionary distance between two
#' genomes.  The larger the value, the more divergent the genomes.
#'
#' @family Nucleotide distances
#'
#' @author Charles Plessy
#'
#' @references May, A. C. W. (2004). Percent Sequence Identity: The Need to Be Explicit.
#' *Structure*, 12(5), 737–738. \doi{10.1016/j.str.2004.04.001}
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' P_distance(parameters)
#'
#' @export

P_distance <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  m <- m[c("A", "C", "G", "T"), c("A", "C", "G", "T")]
  P <- prop.table(m)
  1 - sum(diag(P))
}
