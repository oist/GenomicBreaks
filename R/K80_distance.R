#' Kimura's 2-parameter distance
#'
#' The Kimura's 2-parameter (K80) distance distinguishes between transitions and transversions.  It is calculated as \eqn{-\frac{1}{2} \ln\left(\left(1 - 2p - q\right) \times \sqrt{1 - 2q}\right)}, where \eqn{p} is the probability of transition and \eqn{q} is the probability of transversion. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#K80_model_(Kimura_1980)) for more details.
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @family Nucleotide distances
#'
#' @author Zikun Yang
#'
#' @references Kimura, M. (1980). "A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences." *Journal of Molecular Evolution*, 16, 111–120. DOI: [10.1007/BF01731581](https://doi.org/10.1007/BF01731581)
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' K80_distance(exampleSubstitutionMatrix)
#'
#' @export

K80_distance <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  m <- m[c("A", "C", "G", "T"), c("A", "C", "G", "T")]
  P <- prop.table(m)
  p <- P["A", "G"] + P["G", "A"] +
       P["C", "T"] + P["T", "C"]
  q <- P["A", "C"] + P["C", "A"] +
       P["G", "T"] + P["T", "G"] +
       P["A", "T"] + P["T", "A"] +
       P["C", "G"] + P["G", "C"]
  - 0.5 * log((1 - 2 * p - q) * sqrt(1 - 2 * q))
}
