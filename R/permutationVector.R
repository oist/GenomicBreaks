#' Permutation vector
#'
#' Represent _GenomicBreaks_ objects as a _permutation vector_.
#'
#' In a `GBreaks` object properly sorted by its _target_ ranges, the _query_
#' ranges …
#'
#' @param gb A [`GBreaks`] object.
#'
#' @returns A vector of non-zero integers representing the position of the
#' genomic content in the `GBreaks` object.
#'
#' @export

permutationVector <- function(gb) {
  # Transform query ranges into simple indices
  o <- order(order(gb$query))
  o * ifelse(strand(gb) == '-', -1L, 1L)
}
