#' Flag translocations
#'
#' Flag ranges that start a triplet that would be colinear if the central pair
#' were colinear, and that are not an inversion.
#'
#' The result of the flagging differs if the _target_ and _query_ ranges are
#' swapped.  Therefore, by default this function will search for translocations
#' on both cases, and flag both ranges.
#'
#' Here is a trivial example of a translocation:
#'
#' ```
#' ┌──────────────┬──────────────┬──────────────┐
#' │ chrA:101-200 │ chrA:201-300 │ chrA:301-400 │ (Target genome)
#' └──────────────┴──────────────┴──────────────┘
#'        +               +             +         (Alignment direction)
#' ┌──────────────┬──────────────┬──────────────┐
#' │ chrB:101-200 │ chrC:201-300 │ chrB:301-400 │ (Query genome)
#' └──────────────┴──────────────┴──────────────┘
#' ```
#'
#' @param gb A `GBreaks()` object.
#' @param tol Tolerance window for the distance between two ranges.
#' @param both Flag both the _target_ and _query_ ranges.
#'
#' @return Returns the `GBreaks` object with an extra `tra` metadata column.
#'
#' @family Flagging functions
#' @family Translocation functions
#' @family Structural variants
#'
#' @examples
#' flagTranslocations(exampleTranslocation)
#' flagTranslocations(exampleTranslocation2)
#' plotApairOfChrs(exampleTranslocation)
#' flagTranslocations(exampleDeletion)
#' flagTranslocations(exampleInsertion)
#' flagTranslocations(exampleInsertion  |> swap(sort = TRUE))
#' flagTranslocations(exampleInsertion) |> swap(sort = TRUE)
#' flagTranslocations(sort(reverse(exampleDeletion)))
#' flagTranslocations(exampleInversion)
#' flagTranslocations(exampleColinear3)
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede
#' @export

flagTranslocations <- function (gb, tol = Inf, both = TRUE) {
  if(length(gb) == 0) {
    gb$tra <- Rle(logical(0))
    return(gb)
  }
  if(length(gb) == 1) {
    gb$tra <- Rle(FALSE)
    return(gb)
  }
  # Enforce sorting, to guarantee colinearity on the _target_ranges
  if (isFALSE(isSorted(gb))) stop ("Can not run on non-sorted objects.")
  gb.bak <- gb # save the original object
  # Calculate distance to next and next-next entry
  # (here, "next" is next in order; do not confuse with preceding/following)
  gb <- flagColinearAlignments(gb, details = TRUE) # Provides qfoll and qprev for later
  gb <- dist2next(gb, step = 2, ignore.strand = TRUE)
  gb$tdist2 <- gb$tdist
  gb$qdist2 <- gb$qdist
  gb <- dist2next(gb, step = 1, ignore.strand = TRUE) # Overrides strand-aware values from flagColinearAlignments
  gb$twoSameStrand <- strand(gb) == c(strand(tail(gb, -2)), factor(c("*", "*")))

  gb.bak$tra <-
    # There is a distance, therefore they are on the same sequence.
    (! (is.na(gb$tdist) | is.na(gb$tdist2)))  &
    # Query colinear with the next-next entry…
    (
      ( decode(strand(gb)) %in% c("+", "*") & gb$qfoll == 2 ) |
      ( decode(strand(gb)) %in%    c("-") & gb$qprev == 2 )
    ) &
    # … colinear if same strand
    gb$twoSameStrand

  if(tol < Inf) stop("Not implemented yet")

  gb.bak$tra[is.na(gb.bak$tra)] <- FALSE

  if (isTRUE(both)) {
    gb.swapped <- swap(gb.bak, sort = TRUE)
    gb.swapped <- flagTranslocations(gb.swapped, tol = tol, both = FALSE)
    gb.swapped$query$tra <- gb.swapped$tra # Save the flag on the _query_ because swap will drop the _target_ one.
    gb.reswapped <- swap(gb.swapped, sort = TRUE)
    gb.bak$query$tra <- gb.reswapped$tra
  }

  gb.bak
}

#' Show translocations and their flanking blocks.
#'
#' @param gb A [`GBreaks`] object processed with [`flagTranslocations`].
#'
#' @param rename Replace range names by their numeric order before subsetting.
#'
#' @family Translocation functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' part of a translocation triplet have been discarded.  If the object was missing
#' the `tra` metadata column, return the object after discarding all of its
#' ranges.
#'
#' @examples
#' showTranslocations(flagTranslocations(exampleTranslocation))
#'
#' @export

showTranslocations <- function(gb, rename = TRUE) {
  if (is.null(gb$tra)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  flagPos <- which(gb$tra) + 1
  flagContext <- c(flagPos -1 , flagPos, flagPos + 1) |> unique() |> sort()
  gb[flagContext]
}

#' Extract central blocks in translocations
#'
#' @param gb A [`GBreaks`] object processed with [`flagTranslocations`].
#'
#' @param rename Replace range names by their numeric order before subsetting.
#'
#' @param remove Filter out instead of filtering in.
#'
#' @family Translocation functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' the central part of an inversion triplet have been discarded.  If the object
#' was missing the `tra` metadata column, return the object after discarding all
#' of its ranges.
#'
#' @examples
#' filterTranslocations(flagTranslocations(exampleTranslocation))
#' filterTranslocations(remove = TRUE, flagTranslocations(exampleTranslocation))
#'
#' @export

filterTranslocations <- function(gb, rename = TRUE, remove = FALSE) {
  if (is.null(gb$tra)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  flagPos <- which(gb$tra) + 1
  if (isTRUE(remove))
    flagPos <- -(flagPos)
  gb[flagPos]
}
