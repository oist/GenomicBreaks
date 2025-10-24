#' GC pressure
#'
#' Computes the following equation:
#' \deqn{\frac{\text{W}\to\text{S} \quad-\quad \text{S}\to\text{W}}{\text{W}\to\text{S} \quad+\quad \text{S}\to\text{W}}}{GC* = (W→S - S→W) / (W→S + S→W)}
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @note This index was suggested to by by M365 Copilot while using it to write
#' the [`GCequilibrium`] function and its documentation.  Despite my attempts I
#' could not find a paper that uses it explicitely…  I am including it in this
#' package because I am screening as many indices that I can extract from an
#' alignment pair.
#'
#' @returns A numeric value between \eqn{-1} and \eqn{1}.
#'
#' @family Alignment statistics
#'
#' @author Charles Plessy
#'
#' @examples
#' GCpressure(exampleSubstitutionMatrix)
#' GCpressure(t(exampleSubstitutionMatrix))
#'
#' @export

GCpressure <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  W <- c("A", "T")
  S <- c("G", "C")
  a <- sum(m[W, S])
  b <- sum(m[S, W])
  (a - b) / (a + b)
}
