#' The GenomicBreaks class
#'
#' The `GBreaks` class is a simple wrapper to the `GRanges` class
#'
#' @export

GBreaks <- setClass("GBreaks", contains = "GRanges")


#' Conversion from [`CNEr::Axt`] objects
#'
#' @importFrom CNEr first second
#' @name as
#' @export

setAs("Axt", "GBreaks", function(from) {

  # First genome of axt object is target genome, second is query.
  gb       <- granges(first(from))
  gb$query <- granges(second(from))

  # In GBreaks object, strand information is carried by the target genome ranges
  strand(gb) <- strand(gb$query)
  strand(gb$query) <- "+"

  GBreaks(gb)
})
