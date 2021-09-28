#' Compute all flags
#'
#' Flag ranges that are at the beginning of a colinear duple, or an inversion,
#' insertion or interruption triplet.
#'
#' @param gb A `GBreaks` object.
#'
#' @return Returns the `GBreaks` object with an extra `flag` metadata column.
#'
#' @family Flagging functions
#'
#' @examples
#' flagAll(exampleInversion)
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
