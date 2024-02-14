#' Breakpoint graph
#'
#' Represent _GenomicBreaks_ objects as _breakpoint graphs_.
#'
#' In a `GBreaks` object properly sorted by its _target_ ranges, the _query_
#' ranges …
#'
#' @param gb A [`GBreaks`] object.
#'
#' @returns A vector of positive integers representing the position of the
#' breakpoints in the `GBreaks` object.
#'
#' @export

bpGraph <- function(gb) {
  # Transform query ranges into simple indices
  o <- order(order(gb$query))
  o <- o * ifelse(strand(gb) == '-', -1L, 1L)
  bpg <- sapply(o, function(i) {
    if (i>0) { c( 2L * i,        2L * i + 1L) }
    else     { c(-2L * i + 1L , -2L * i     ) }
  })
  c(1L, bpg, 2L * length(o) + 2L)
}
