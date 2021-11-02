#' Compute all flags
#'
#' Flag ranges that are at the beginning of a colinear duplet, or an inversion
#' or a translocation triplet.
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
#' @include flagColinearAlignments.R Translocations.R flagInversions.R
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

  Tra <- flagTranslocations(gb)$tra |> decode()
  stopifnot(all(is.na(gb$flag[Tra])))
  gb$flag[Tra] <- "Tra"

  Dbl <- flagDoubleInversions(gb)$Dbl |> decode()
  stopifnot(all(is.na(gb$flag[Dbl])))
  gb$flag[Dbl] <- "Dbl"

  gb
}
