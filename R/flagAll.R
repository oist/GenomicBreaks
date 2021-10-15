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

  Col <- flagColinearAlignments(gb)$colinear |> decode()
  stopifnot(all(is.na(gb$flag[Col])))
  gb$flag[Col] <- "Col"

  Inv <- flagInversions(gb)$inv |> decode()
  stopifnot(all(is.na(gb$flag[Inv])))
  gb$flag[Inv] <- "Inv"

  Ins <- flagInsersions(gb)$ins |> decode()
  stopifnot(all(is.na(gb$flag[Ins])))
  gb$flag[Ins] <- "Ins"

  Tra <- flagTranslocations(gb)$tra |> decode()
  stopifnot(all(is.na(gb$flag[Tra])))
  gb$flag[Tra] <- "Tra"

  gb
}
