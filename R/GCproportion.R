#' Proportion of GC in alignment
#'
#' The proportion of GC in the _target_ genome in an alignment.  For the same
#' on the _query_ genome, just transpose the input matrix.
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @returns A numeric value between 0 and 1.
#'
#' @family Nucleotide distances
#'
#' @author Charles Plessy
#'
#' @examples
#' GCproportion(exampleSubstitutionMatrix)
#' GCproportion(t(exampleSubstitutionMatrix))
#'
#' @export

GCproportion <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  sum(colSums(m)[c("G","C")]) / sum(colSums(m)[c("A", "C", "G","T")])
}
