#' Swap target and query genomes
#'
#' Produce a new object in which the information about the _target_ and
#' _query_ genomes have been inverted.  As a large number of _GenomicBreaks_
#' functions require `GBreaks` objects that were sorted with strand information
#' ignored, a `sort` option is provided to do that.
#'
#' @param gb A [`GBreaks`] object.
#' @param sort Sort the returned object, ignoring strand.
#'
#' @return A new `GBreaks` object representing the input's _query_ ranges,
#' with the input's _target_ ranges in the _query_ slot.  The strand information
#' is unchanged.
#'
#' @family modifier functions
#'
#' @examples
#' swap(exampleDeletion)
#' swap(exampleDeletion, sort = TRUE)
#'
#' @export

swap <- function(gb, sort = FALSE) {
  gb2               <- gb$query
  gb2$query         <- granges(gb)
  strand(gb2)       <- strand(gb)
  strand(gb2$query) <- "*"
  if (isTRUE(sort))
    gb2 <- sort(gb2, ignore.strand = TRUE)
  GBreaks(gb2)
}
