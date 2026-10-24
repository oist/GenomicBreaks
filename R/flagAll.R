#' Compute all flags
#'
#' Flag ranges that are at the beginning of a colinear duplet, or an inversion
#' or a translocation triplet.
#'
#' @param gb A `GBreaks` object.
#'
#' @section Unstranded ranges (experimental):
#' On ranges of the `*` strand, a range can be both colinear with the next
#' one and the start of a payload region of [`flagTranslocations()`]; it is
#' then flagged as colinear.  This behaviour may change in future versions.
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
  if (length(gb) == 0) {
    gb$flag <- logical(0)
    return(gb)
  }
  gb$flag <- NA

  Col <- flagColinearAlignments(gb)$colinear |> decode()
  stopifnot(all(is.na(gb$flag[Col])))
  gb$flag[Col] <- "Col"

  Inv <- flagInversions(gb)$inv |> decode()
  stopifnot(all(is.na(gb$flag[Inv])))
  gb$flag[Inv] <- "Inv"

  Tra <- flagTranslocations(gb)$tra |> decode()
  # Unstranded ranges (experimental): adjacency is tested in both directions,
  # so that a range can be colinear with the next one and also start a
  # payload region.  Colinearity takes precedence.
  Tra[Tra & decode(strand(gb)) == "*" & !is.na(gb$flag)] <- FALSE
  stopifnot(all(is.na(gb$flag[Tra])))
  gb$flag[Tra] <- "Tra"

  Dbl <- flagDoubleInversions(gb)$Dbl |> decode()
  stopifnot(all(is.na(gb$flag[Dbl])))
  gb$flag[Dbl] <- "Dbl"

  Twi <- flagTwinInversions(gb)$twi |> decode()
  stopifnot(all(is.na(gb$flag[Twi])))
  gb$flag[Twi] <- "Twi"

  gb
}
