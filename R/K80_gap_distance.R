#' Kimura's 2-parameter distance extended with gaps
#'
#' Nishimaki and Sato's K2P + Gap distance extends Kimura's 2-parameter
#' distance by treating gaps as insertion/deletion states instead of removing
#' gap-containing sites.
#'
#' The distance is calculated as
#' \deqn{
#' K = \frac{3}{4} w \log(w)
#'     - \frac{w}{2} \log\left((S - P) \sqrt{S + P - Q}\right)
#' }
#' where \eqn{S} is the probability of identical nucleotide pairs,
#' \eqn{P} is the probability of transition-type nucleotide pairs,
#' \eqn{Q} is the probability of transversion-type nucleotide pairs, and
#' \eqn{w} is the nucleotide occupancy probability.
#'
#' When there are no gaps, \eqn{w = 1}, and this expression reduces to the
#' usual K80 distance.
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome. As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @author Charles Plessy
#'
#' @references
#' Kimura, M. (1980). "A simple method for estimating evolutionary rates of base
#' substitutions through comparative studies of nucleotide sequences."
#' *Journal of Molecular Evolution*, 16, 111–120.
#' DOI: \doi{10.1007/BF01731581}
#'
#' Nishimaki, T. and Sato, K. (2019). "An Extension of the Kimura
#' Two-Parameter Model to the Natural Evolutionary Process."
#' *Journal of Molecular Evolution*, 87, 60–67.
#' DOI: \doi{10.1007/s00239-018-9885-1}
#'
#' @returns Returns a numeric value showing the evolutionary distance between
#' two genomes. The larger the value, the more different the two genomes are.
#'
#' @examples
#' K80_gap_distance(exampleSubstitutionMatrix)
#'
#' # When there are no gaps, it returns the same as the K80 distance
#' nogaps <- exampleSubstitutionMatrix
#' nogaps["-",] <- 0
#' nogaps[,"-"] <- 0
#' K80_gap_distance(nogaps)
#' K80_distance(nogaps)
#'
#' @export

K80_gap_distance <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  nuc <- c("A", "C", "G", "T")
  states <- c(nuc, "-")
  m <- m[states, states, drop = FALSE]
  if (all(m == 0)) return(NA)
  P <- prop.table(m)

# identical nucleotide pairs.
  s <- sum(diag(P[nuc, nuc, drop = FALSE]))

# transition-type nucleotide pairs.
  p <- P["A", "G"] + P["G", "A"] +
       P["C", "T"] + P["T", "C"]

# transversion-type nucleotide pairs.
  q <- P["A", "C"] + P["C", "A"] +
       P["G", "T"] + P["T", "G"] +
       P["A", "T"] + P["T", "A"] +
       P["C", "G"] + P["G", "C"]

# nucleotide-gap pairs.
  g <- sum(P[nuc, "-", drop = TRUE]) +
       sum(P["-", nuc, drop = TRUE])

# nucleotide occupancy probability (a nucleotide-gap site contributes 1/2)
  w <- s + p + q + g / 2
  if (w <= 0) return(NA_real_)

  # As with standard K80, the logarithm becomes undefined for saturated or
  # incompatible count patterns.
  a <- s - p
  b <- s + p - q
  if (a <= 0 || b <= 0) return(NaN)

  0.75 * w * log(w) - 0.5 * w * log(a * sqrt(b))
}
