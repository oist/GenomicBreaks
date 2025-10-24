#' Proportion of gaps in alignment
#'
#' The proportion of nucleotides that are part of an alignment but not matched
#' to another nucleotide.  If you divide it by the percent difference, you obtain
#' the S/I ratio studied by Chen and coll., 2009.
#'
#' @param m A matrix of counts or probabilities for bases of the _target_ genome
#' to be aligned to bases on the _query_ genome.  As a convenience it can also
#' receive a list produced by the [`readTrainFile()`] function, containing this
#' matrix.
#'
#' @note This is not a measure of distance. At greater evolutionary divergence,
#' genome alignment becomes increasingly difficult, leaving only the most
#' conserved regions alignable. Consequently, the proportion of gaps tends to
#' decrease.
#'
#' @returns A numeric value between 0 and 1.
#'
#' @family Nucleotide distances
#'
#' @author Charles Plessy
#'
#' @references Chen JQ, Wu Y, Yang H, Bergelson J, Kreitman M, Tian D. (2009).
#' Variation in the ratio of nucleotide substitution and indel rates across
#' genomes in mammals and bacteria. *Mol Biol Evol.* 26(7):1523-31.
#' \doi{10.1093/molbev/msp063}
#'
#' @examples
#' parameters <- readTrainFile(system.file("extdata/example.train", package = "GenomicBreaks"))
#' gapProportion(parameters)
#'
#' @export

gapProportion <- function(m) {
  if(is.list(m)) m <- m$probability_matrix
  if (! '-' %in% c(colnames(m), rownames(m))) return(0)
  (sum(m['-',]) + sum(m[,'-'])) / sum(m)
}
