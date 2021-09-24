#' Compute all flags
#'
#' Flag ranges that are at the beginning of a colinear duplet, or an inversion,
#' insertion or interruption triplet.
#'
#' @param gb A `GenomicBreaks` object.
#'
#' @return Returns the `GenomicBreaks` object with an extra `flag` metadata
#' column.
#'
#' @family Flagging functions
#'
#' @examples
#' inv       <- GRanges(c("XSR:101-200:+", "XSR:201-300:-",  "XSR:301-400:+"))
#' inv$query <- GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400"))
#' flagAll(inv)
#'
#' @include flagColinearAlignments.R flagInsersions.R flagInversions.R
#' @importFrom S4Vectors decode
#' @export

flagAll <- function (gb) {
  gb$flag <- NA
  gb$flag[flagColinearAlignments(gb)$colinear |> decode()] <- "Col"
  gb$flag[flagInversions(gb)$inv              |> decode()] <- "Inv"
  gb$flag[flagInsersions(gb)$ins              |> decode()] <- "Ins"
  gb
}
