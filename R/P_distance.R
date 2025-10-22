#' Percent difference Distance
#'
#' The P-distance is simply the percentage of the aligned base pairs where the
#' two bases are different.  See Details for a discussion on the many ways to
#' compute percent identity (and therefore difference.)
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#' @param denominator Denominator according to the nomenclature of May (2004).
#' Default is `L3`.
#'
#' @details The P-distance is computed as 1 minus the percent sequence identity
#' according to the definitions of May (2004), which is the ratio between the
#' number of matching bases and four possible denominators in use at that time.
#'
#'  - L1: Length of the shorter sequence.
#'  - L2: Number of aligned positions, i.e., alignment length (includes gaps, if any).
#'  - L3: Number of aligned residue pairs, i.e., identities and nonidentities (excludes gaps, if any).
#'  - L4: Arithmetic mean sequence length.
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

P_distance <- function(m, denominator = c("L3", "L1", "L2", "L4")) {
  denominator <- match.arg(denominator)
  if(is.list(m)) m <- m$probability_matrix
  non_gap <- c("A", "C", "G", "T")
  numerator   <- sum(diag(m[non_gap, non_gap]))
  denominator <- switch( denominator
                       , L1 =  min(sum(m[non_gap,        ]),
                                   sum(m[       , non_gap]))
                       , L2 =      sum(m)
                       , L3 =      sum(m[non_gap, non_gap])
                       , L4 = mean(sum(m[non_gap,        ]),
                                   sum(m[       , non_gap]))
                       )
  1 - numerator / denominator
}
