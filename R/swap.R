#' Swap target and query genomes
#'
#' Produce a new object in which the information about the _target_ and
#' _query_ genomes have been inverted.
#'
#' @param gb A [`GBreaks`] object
#'
#' @return A new `GBreaks` object representing the input's _query_ ranges,
#' with the input's _target_ ranges in the _query_ slot.  The strand information
#' is unchanged.
#'
#' @family modifier functions
#'
#' @examples
#' swap(exampleColinear3)
#'
#' @export

swap <- function(gb) {
  gb2               <- gb$query
  gb2$query         <- granges(gb)
  strand(gb2)       <- strand(gb)
  strand(gb2$query) <- "*"
  GBreaks(gb2)
}
