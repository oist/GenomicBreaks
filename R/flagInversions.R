#' Flag inversions
#'
#' Flag ranges that start a triplet that would be colinear if the central pair
#' were aligned to the opposite strand.
#'
#' Here is a trivial example of an inversion.
#'
#' ```
#' ┌──────────────┬──────────────┬──────────────┐
#' │ chrA:101-200 │ chrA:201-300 │ chrA:301-400 │ (Target genome)
#' └──────────────┴──────────────┴──────────────┘
#'        +               -             +         (Alignment direction)
#' ┌──────────────┬──────────────┬──────────────┐
#' │ chrB:101-200 │ chrB:201-300 │ chrB:301-400 │ (Query genome)
#' └──────────────┴──────────────┴──────────────┘
#' ```
#'
#' @param gb A `GBreaks()` object.
#' @param tol Tolerance window for the distance between two ranges.
#'
#' @return Returns the `GBreaks` object with an extra `inv` metadata column.
#'
#' @family Flagging functions
#' @family Inversion functions
#' @family Structural variants
#'
#' @examples
#' flagInversions(exampleInversion)
#' flagInversions(exampleInversion, tol = 19)
#' plotApairOfChrs(exampleInversion)
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede
#' @export

flagInversions <- function (gb, tol = Inf) {
  if (isFALSE(isSorted(gb))) stop ("Can not run on non-sorted objects.")
  gb.bak <- gb # save the original object
  gb <- dist2next(gb, ignore.strand = TRUE)
  gb$tdist2 <- c(gb$tdist[-1], Inf)
  gb$qdist2 <- c(gb$qdist[-1], Inf)
  gb$oneNext <- precede(gb$query, tail(gb$query, -1), ignore.strand = T) - seq_along(gb) == 0
  gb$twoNext <- c(gb$oneNext[-1], NA)
  gb$oneDiffStrand <- strand(gb) != c(strand(tail(gb, -1)), factor("*"))
  gb$twoSameStrand <- strand(gb) == c(strand(tail(gb, -2)), factor(c("*", "*")))
  gb.bak$inv <- gb$oneNext          &
            gb$twoNext          &
            gb$oneDiffStrand    &
            gb$twoSameStrand    &
            gb$tdist  < tol + 1 &   # We need `< tol + 1`
            gb$tdist2 < tol + 1 &   # and not `<= tol`
            gb$qdist  < tol + 1 &   # Because of comparisons to
            gb$qdist2 < tol + 1     # infinite distance returned by dist2next.
  if (any(is.na(gb.bak$inv)))
    gb.bak[is.na(gb.bak$inv)]$inv <- FALSE
  gb.bak
}

#' Show inversions and their flanking blocks.
#'
#' @param gb A [`GBreaks`] object processed with [`flagInversions`].
#'
#' @param rename Replace range names by their numeric order before subsetting.
#'
#' @family Inversion functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' part of an inversion triplet have been discarded.  If the object was missing
#' the `inv` metadata column, return the object after discarding all of its
#' ranges.
#'
#' @examples
#' showInversions(flagInversions(exampleInversion))
#'
#' @export

showInversions <- function(gb, rename = TRUE) {
  if (is.null(gb$inv)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  invPos <- which(gb$inv) + 1
  invContext <- c(invPos -1 , invPos, invPos + 1) |> unique() |> sort()
  gb[invContext]
}

#' Extract central blocks in inversions
#'
#' @param gb A [`GBreaks`] object processed with [`flagInversions`].
#'
#' @param rename Replace range names by their numeric order before subsetting.
#'
#' @family Inversion functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' the central part of an inversion triplet have been discarded.  If the object
#' was missing the `inv` metadata column, return the object after discarding all
#' of its ranges.
#'
#' @examples
#' filterInversions(flagInversions(exampleInversion))
#'
#' @export

filterInversions <- function(gb, rename = TRUE) {
  if (is.null(gb$inv)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  invPos <- which(gb$inv) + 1
  invContext <- c(invPos) |> unique() |> sort()
  gb[invContext]
}

#' Isolate the left-side gaps in inversions
#'
#' @param gb A [`GBreaks`] object.
#'
#' @return Returns a [`GRanges`] object representing the left-side gaps in the
#' `GBreaks` object.
#'
#' @family Inversion functions
#'
#' @examples
#' leftInversionGaps(flagInversions(exampleInversion))
#'
#' @export

leftInversionGaps <- function(gb) {
  # Flag inversions.
  #  WHAT FOLLOWS ASSUMES THAT THE FLAG IS ON THE LEFT-SIDE BLOCK OF THE TRIPLE
  lgaps <- flagInversions(gb)
  # Extract inversions
  invs <- lgaps[lgaps$inv]
  # Extend inversion of 1 bp so that they overlap with their neighbor gap
  invs <- shift(invs, 1)
  # Remove strand information
  strand(lgaps) <- "*"
  # Then extract gap positions ignoring strand
  lgaps <- gaps(lgaps, start=min(start(lgaps)), end=max(end(lgaps)))
  # Remove bogus gaps on + and - strands
  lgaps <- lgaps[strand(lgaps) == "*"]
  # Return the gaps overlapping with the flagged inversions
  subsetByOverlaps(lgaps, invs)
}
