#' Jukes-Cantor 1969 distance
#'
#' The Jukes and Cantor 1969 (JC69) distance is calculated as \eqn{-\frac{3}{4} \ln(1 - \frac{4}{3} p )}, where \eqn{p} is the sum of the probabilities of the alignment. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) for more details.
#'
#' @references Jukes, T.H. & Cantor, C.R. (1969). "Evolution of protein molecules." In *Mammalian Protein Metabolism* (pp. 21–132). Academic Press.
#'
#' @param m Two possible types of input are available:
#' 1. A matrix of counts or probabilities for bases of the _target_ genome to be aligned to bases on the _query_ genome;
#' 2. The mismatch count (a positive integer) or proportion (a number between `0.0` and `1.0`).
#' @param tot The total count of nucleotides, in case the parameter `m` is the _count_ of mismatches.
#' @param adjust_p A boolean flag. If `TRUE`, the distance is scaled between `0` and `0.75` to ensure the logarithm stays positive.
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @author Zikun Yang
#' @author Priscila Biller
#'
#' @returns Returns a numeric value representing the evolutionary distance between two genomes. The greater the value, the more genetically different the genomes are.
#'
#' @examples
#' # Using a substitution matrix as input:
#' d <- JC69_distance(exampleSubstitutionMatrix)
#'
#' # Using numerical values as input:
#' d <- JC69_distance(25, 100)
#'
#' @export
JC69_distance <- function(m, tot=1.0, adjust_p=FALSE) {
  UseMethod("JC69_distance")
} 

#' @export
JC69_distance.numeric <- function(m, tot=1.0, adjust_p=FALSE) {
  # The p-distance, i.e., the fraction of nucleotides that are different.
  p <- m / tot
  # In this model, at equilibrium, the p-distance between two sequences
  # that began as identical should be 75%, and should never exceed 75%.
  # This 'trick' ensures the log argument stays valid (and the expression
  # remains well-behaved) even if the raw distance would be above 75%.
  if (adjust_p == TRUE) { p <- p * 0.75 }
  # Case where log is negative: the distance is undefined.
  else if (p > 0.75) {return(NA)}
  - 0.75 * log(1 - 4 * p / 3)
}

#' @export
JC69_distance.matrix <- function(m, ...) {
  p <- P_distance(m, denominator = "L3")
  JC69_distance(p, ...)
}
