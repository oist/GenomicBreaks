#' Felsenstein's 1981 Distance
#'
#' The Felsenstein 1981 (F81) distance corrects for multiple substitutions like
#' the [`JC69_distance()`] function, but also accounts for GC content. It is
#' calculated as \eqn{-E \ln\left(\frac{E - p}{E}\right)}, where \eqn{E} equals
#' \eqn{1} minus the sum of the squares of the average nucleotide frequencies
#' (McGuire and coll., 1999).
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
#' @references
#' Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum
#' likelihood approach.  *Journal of Molecular Evolution*, 17, 368–376.
#' \doi{10.1007/BF01734359}
#'
#' McGuire, G., Prentice, M. J., & Wright, F. (1999). Improved error bounds for
#' genetic distances from DNA sequences.  *Biometrics*, 55(4), 1064–1070.
#' \doi{10.1111/j.0006-341x.1999.01064.x}
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' F81_distance(parameters)
#'
#' @export

F81_distance <- function(train_parameters) {
  P <- train_parameters$probability_matrix
  target_freqs <- rowSums(P)
  query_freqs  <- colSums(P)
  avg_freqs    <- (target_freqs + query_freqs) / 2
  p <- 1 - sum(diag(P))
  E <- 1 - sum(avg_freqs^2)
  if (p >= E) return(Inf)
  -E * log((E - p) / E)
}
