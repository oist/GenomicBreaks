#' Jukes-Cantor 1969 distance
#'
#' The Jukes and Cantor 1969 (JC69) distance is calculated as \eqn{-\frac{3}{4} \ln(1 - \frac{4}{3} p )}, where \eqn{p} is the sum of the probabilities of the alignment. See the [Wikipedia](https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_(Jukes_and_Cantor_1969)) for more details.
#'
#' @references Jukes, T.H. & Cantor, C.R. (1969). "Evolution of protein molecules." In *Mammalian Protein Metabolism* (pp. 21–132). Academic Press.
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @family Alignment statistics
#' @family Similarity indexes
#'
#' @author Zikun Yang
#'
#' @returns Returns a numeric value show the evolutionary distance between two genomes. the larger the value, the more different the two genomes are.
#'
#' @examples
#' JC69_distance(exampleSubstitutionMatrix)
#'
#' @export

JC69_distance <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  m <- m[c("A", "C", "G", "T"), c("A", "C", "G", "T")]
  if (all(m == 0)) return(NA)
  P <- prop.table(m)
  p <- 1 - sum(diag(P))
  - 0.75 * log(1 - 4 * p / 3)
}
