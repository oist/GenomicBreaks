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
#' flagInversions(exampleInversion |> reverse() |> sort(ignore.strand = TRUE))
#' plotApairOfChrs(exampleInversion |> reverse())
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede follow
#' @export

flagInversions <- function (gb, tol = Inf) {
  if (isFALSE(isSorted(gb))) stop ("Can not run on non-sorted objects.")
  gb.bak <- gb # save the original object
  gb <- dist2next(gb, ignore.strand = TRUE)
  gb$tdist2 <- c(gb$tdist[-1], Inf)
  gb$qdist2 <- c(gb$qdist[-1], Inf)
  # Check that the next and next-next are colinear, ignoring strand
  gb$oneNext <- precede(gb$query, tail(gb$query, -1), ignore.strand = T) - seq_along(gb) == 0
  gb$twoNext <- c(gb$oneNext[-1], NA)
  gb$colPlus <- strand(gb) %in% c("+", "*") &  gb$oneNext & gb$twoNext
  # Check that the next and next-next are anti-colinear, ignoring strand
  gb$onePrev <- follow(gb$query, tail(gb$query, -1), ignore.strand = T) - seq_along(gb) == 0
  gb$twoPrev <- c(gb$onePrev[-1], NA)
  gb$colMin  <- strand(gb) %in%     "-"     &  gb$onePrev & gb$twoPrev
  # Is the current strand the same as the next ?
  gb$nextDiffStrand     <- strand(gb) != c(strand(tail(gb, -1)), factor(    "*"    ))
  # Is the current strand the same as the next-next ?
  gb$nextNextSameStrand <- strand(gb) == c(strand(tail(gb, -2)), factor(c("*", "*")))
  gb.bak$inv <-
           (gb$colPlus | gb$colMin)  &
            gb$nextDiffStrand        &
            gb$nextNextSameStrand    &
            gb$tdist  < tol + 1 &   # We need `< tol + 1`
            gb$tdist2 < tol + 1 &   # and not `<= tol`
            gb$qdist  < tol + 1 &   # Because of comparisons to
            gb$qdist2 < tol + 1     # infinite distance returned by dist2next.
  if (any(is.na(gb.bak$inv)))
    gb.bak[is.na(gb.bak$inv)]$inv <- FALSE
  gb.bak
}

#' Flag double inversions
#'
#' Two consecutive and overlapping inversions will generate patterns that can
#' this function aims to detect.
#'
#' @param gb A [`GBreaks()`] object.
#' @param details Returns more metadata columns if `TRUE`.
#'
#' @family Flagging functions
#' @family Inversion functions
#' @family Structural variants
#'
#' @examples
#' # Start colinear.  Lower case meands minus strand
#' z0 <- GBreaks(target = GRanges(c(A="T:10-15:+", B="T:20-25:+", C="T:30-35:+")),
#'               query  = GRanges(c(A="Q:10-15",   B="Q:20-25",   C="Q:30-35")))
#'
#' # Swap coordinates of A and B on the query and flip strands  ABC -> baC
#' z1 <- GBreaks(target = GRanges(c(A="T:10-15:-", B="T:20-25:-", C="T:30-35:+")),
#'               query  = GRanges(c(a="Q:20-25",   b="Q:10-15",   C="Q:30-35")))
#'
#' # Now query order is b - a - C.  Swap a and C and flip strands baC -> bcA
#' z2 <- GBreaks(target = GRanges(c(A="T:10-15:+", B="T:20-25:-", C="T:30-35:-")),
#'               query  = GRanges(c(A="Q:30-35",   b="Q:10-15",   c="Q:20-25")))
#'
#' # Altogether, there are:
#' # ABC -> baC -> bcA
#' # ABC -> Acb -> Cab
#' # cba -> BCa -> BAc
#' # cba -> cAB -> aCB
#'
#'
#' # z2 has same topoloty as package example
#' plotApairOfChrs(z2)
#' plotApairOfChrs(exampleDoubleInversion1)
#'
#' flagDoubleInversions(exampleDoubleInversion1 )
#'
#' exampleDoubleInversion1Rev <- reverse(exampleDoubleInversion1) |> sort(ignore.strand = TRUE)
#' exampleDoubleInversion1Rev[2:4] |> plotApairOfChrs()
#'
#' @author Charles Plessy
#'
#' @export

flagDoubleInversions <- function(gb, details = FALSE) {
  if (isFALSE(isSorted(gb))) stop ("Can not run on non-sorted objects.")

  gb.bak <- gb # Backup the object

  # Look ahead for the two next alignments
  gb$qnext0  <- precede(gb$query) - seq_along(gb$query)
  gb$qnext1  <- c(tail(gb$qnext0, -1), NA)
  gb$qnext2  <- c(tail(gb$qnext0, -2), NA, NA)

  gb$strand0 <- strand(gb)
  gb$strand1 <- c(tail(strand(gb), -1), factor("*"))
  gb$strand2 <- c(tail(strand(gb), -2), factor("*"), factor("*"))

  # Pattern 1: bcA
  # case (strand)-insensitive order of bcA is 3, 1, 2
  # In that order, strands are +, -, -
  # precede(IPos(c(3,1,2))) - 1:3
  # We do not test value for the 3rd position (NA)
  gb$pat_bcA <-
    ( gb$qnext1  ==  1  &
      gb$qnext2  == -2  &
      gb$strand0 == "+" &
      gb$strand1 == "-" &
      gb$strand2 == "-"   ) |> sapply(isTRUE) # Turn NAs into FALSE

  # Pattern 2: Cab
  # precede(IPos(c(2,3,1))) - 1:3
  gb$pat_Cab <-
    ( gb$qnext0  ==  1  &
      gb$qnext2  == -2  &
      gb$strand0 == "-" &
      gb$strand1 == "-" &
      gb$strand2 == "+"   ) |> sapply(isTRUE)

  # Pattern 3: BAc
  # precede(IPos(c(2,1,3))) - 1:3
  gb$pat_BAc <-
    ( gb$qnext0  ==  2  &
      gb$qnext1  == -1  &
      gb$strand0 == "+" &
      gb$strand1 == "+" &
      gb$strand2 == "-"   ) |> sapply(isTRUE)

  # Pattern 4; aCB
  # precede(IPos(c(1,3,2))) - 1:3
  gb$pat_aCB <-
    ( gb$qnext0  ==  2  &
      gb$qnext2  == -1  &
      gb$strand0 == "-" &
      gb$strand1 == "+" &
      gb$strand2 == "+"   ) |> sapply(isTRUE)

  gb$Dbl <- gb$pat_bcA | gb$pat_Cab |  gb$pat_BAc | gb$pat_aCB

  if (isTRUE(details)) return(gb)

  gb.bak$Dbl <- gb$Dbl

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

#' Extract central blocks in double inversions
#'
#' @param gb A [`GBreaks`] object processed with [`flagDoubleInversions`].
#'
#' @param rename Replace range names by their numeric order before subsetting.
#'
#' @family Inversion functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' the central part of a double inversion quintuplet have been discarded.  If
#' the object was missing the `Dbl` metadata column, return the object after
#' discarding all of its ranges.
#'
#' @examples
#' filterDoubleInversions(flagDoubleInversions(exampleDoubleInversion1))
#'
#' @export

filterDoubleInversions <- function(gb, rename = TRUE) {
  if (is.null(gb$Dbl)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  invPos <- which(gb$Dbl) + 1
  invContext <- c(invPos, invPos + 1, invPos +2) |> unique() |> sort()
  gb[invContext]
}

#' Flip strand of inversions
#'
#' @param gb A [`GBreaks`] object processed with [`flagInversions`].
#'
#' @family Inversion functions
#' @family modifier functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are the
#' central part of an inversion triplet have their strand orientation flipped.
#' The `inv` metadata column is then discarded as it is not valid anymore.  As
#' the former inversion triplets are now collinear, new inversions may possibly
#' found after the object is coalesced again.
#'
#' @examples
#' exampleNestedInversions |> flagInversions()
#' flipInversions(exampleNestedInversions |> flagInversions())
#'
#' @export

flipInversions <- function(gb) {
  if (length(gb) == 0) return(gb)
  if (is.null(gb$inv)) return(gb[0])
  invPos <- which(gb$inv) + 1
  strand(gb)[invPos] <- ifelse(strand(gb)[invPos] == "+", "-", "+")
  gb$inv <- NULL
  gb
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
