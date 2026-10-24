#' Flag translocations
#'
#' Flag the two places where a range that moved between the two genomes is
#' seen in the order of the _target_ genome: the _payload region_, where the
#' moved range is found on the _target_ genome, and the _site_, where it is
#' found on the _query_ genome.
#'
#' A range `A` moved between the two genomes when, for instance:
#'
#' ```
#' Target genome:  X Y . . . U A V
#' Query genome:   X A Y . . . U V
#' ```
#'
#' * The _payload region_ is a triplet of consecutive ranges `U A V` on one
#'   _target_ sequence, such that `U` and `V` are on the same strand and would
#'   be colinear if `A` were removed.  The `tra` flag is set on `U`.
#' * The _site_ is a pair of consecutive ranges `X Y` on one _target_ sequence,
#'   on the same strand, such that `A` is found between them on the _query_
#'   genome, and that `X A Y` would be colinear if `A` were on the same strand
#'   as `X` and `Y`.  The `trs` flag is set on `X`.
#'
#' `A` can be on either strand, and on any _target_ sequence.  When one of the
#' two genomes has no context for `A`, for instance when it is alone on its
#' sequence, only one of the two flags is set.  When both are set, they
#' designate the same range `A`.  The function therefore detects transpositions
#' of single ranges, whether or not they change sequence.
#'
#' Here is a trivial example of a translocation, where the _query_ genome has
#' no context for the moved range, so that only the payload region is flagged:
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
#' @note Before version 0.25.0, the site was flagged only with the `both`
#' option, by a `tra` flag on the _query_ ranges, set on the first range of the
#' site in _query_ order.
#'
#' @param gb A `GBreaks()` object.
#' @param tol Ranges separated by more than `tol` unaligned bases are not
#'        considered adjacent: `U`, `A` and `V` on the _target_ genome and `U`
#'        and `V` on the _query_ genome for payload regions; `X` and `Y` on the
#'        _target_ genome and `X`, `A` and `Y` on the _query_ genome for sites.
#'
#' @return Returns the `GBreaks` object with two extra metadata columns, `tra`
#' and `trs`, flagging the payload regions and the sites.
#'
#' @family Flagging functions
#' @family Translocation functions
#' @family Structural variants
#'
#' @examples
#' flagTranslocations(exampleTranslocation)
#' # The ranges U and V are separated by 100 bases on the query genome.
#' flagTranslocations(exampleTranslocation, tol = 100)
#' flagTranslocations(exampleTranslocation, tol = 99)
#' flagTranslocations(exampleTranslocation2)
#' plotApairOfChrs(exampleTranslocation)
#'
#' # When the genomes are swapped, the moved range has no context on the
#' # target genome, and only the site is flagged.
#' flagTranslocations(swap(exampleTranslocation, sort = TRUE))
#'
#' # A move along one sequence is seen at both places.
#' gb <- permutationVectorToGBreaks(c(1, 3, 4, 2, 5))
#' flagTranslocations(gb)
#'
#' flagTranslocations(exampleDeletion)
#' flagTranslocations(exampleInsertion)
#' flagTranslocations(exampleInsertion |> swap(sort = TRUE))
#' flagTranslocations(sort(reverse(exampleDeletion), ignore.strand = TRUE))
#' flagTranslocations(exampleInversion)
#' flagTranslocations(exampleColinear3)
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede follow
#' @export

flagTranslocations <- function (gb, tol = Inf) {
  if(length(gb) < 2) {
    gb$tra <- gb$trs <- Rle(FALSE, length(gb))
    return(gb)
  }
  # Enforce sorting, to guarantee colinearity on the _target_ ranges
  if (isFALSE(isSorted(gb))) stop ("Can not run on non-sorted objects.")
  gb.bak <- gb # save the original object

  # Offset of the range that each range immediately precedes (qnext0) or
  # follows (qprev0) on the query genome.  The plus-strand tests also apply to
  # unstranded ranges.
  gb$qnext0  <- precede(gb$query) - seq_along(gb$query)
  gb$qprev0  <- follow( gb$query) - seq_along(gb$query)
  plus  <- decode(strand(gb)) %in% c("+", "*")
  minus <- decode(strand(gb)) ==      "-"

  gb$strand0 <- strand(gb)
  gb$strand1 <- c(tail(strand(gb), -1), factor("*"))
  gb$strand2 <- c(tail(strand(gb), -2), factor("*"), factor("*"))

  gb$tseq0   <- seqnames(gb)
  gb$tseq1   <- c(tail(seqnames(gb), -1), factor(NA))
  gb$tseq2   <- c(tail(seqnames(gb), -2), factor(NA), factor(NA))

  # Distances between ranges that must be adjacent (see dist2next()): on the
  # target genome, from range i to i + 1 (tdist0) and from i + 1 to i + 2
  # (tdist1); on the query genome, from range i to range i + 2 (qdist02).
  gb$tdist0  <- dist2next(granges(gb), ignore.strand = TRUE)$dist
  gb$tdist1  <- c(tail(gb$tdist0, -1), Inf)
  gb$qdist02 <- dist2next(granges(gb$query), step = 2)$dist

  # Payload region U A V: U and V are adjacent on the query, that is
  # precede(IPos(c(1, 3, 2))) - 1:3 -> 2, NA, -1 on the plus strand.
  gb$tra <-
    ( gb$tseq2   == gb$tseq0   &
      gb$strand2 == gb$strand0 &
      ( (plus  & gb$qnext0 == 2) |
        (minus & gb$qprev0 == 2) ) &
      gb$tdist0  <= tol &
      gb$tdist1  <= tol &
      gb$qdist02 <= tol ) |> sapply(isTRUE) # Turn NAs into FALSE

  # Site X Y: X precedes the payload A on the query, and A precedes Y, that
  # is precede(IPos(c(1, 3))) - 1:2 -> k, NA with A at position 1 + k, and
  # precede(A) - (1 + k) == 1 - k.  The qnext value of A is looked up through
  # the offset k of X.  Mirror image with follow() on the minus strand.
  gb$qnextPay <- gb$qnext0[seq_along(gb) + gb$qnext0]
  gb$qprevPay <- gb$qprev0[seq_along(gb) + gb$qprev0]
  # Distances on the query genome from X to A and from A to Y.
  pay   <- seq_along(gb) + ifelse(minus, gb$qprev0, gb$qnext0)
  payOK <- ifelse(is.na(pay), seq_along(gb), pay)    # Index without NAs
  next1 <- pmin(seq_along(gb) + 1, length(gb))
  gb$qdistXA <- ifelse(is.na(pay), NA, distance(gb$query,        gb$query[payOK]))
  gb$qdistAY <- ifelse(is.na(pay), NA, distance(gb$query[payOK], gb$query[next1]))
  gb$trs <-
    ( gb$tseq1   == gb$tseq0   &
      gb$strand1 == gb$strand0 &
      ( (plus  & gb$qnext0 != 1 & gb$qnext0 + gb$qnextPay == 1) |
        (minus & gb$qprev0 != 1 & gb$qprev0 + gb$qprevPay == 1) ) &
      gb$tdist0  <= tol &
      gb$qdistXA <= tol &
      gb$qdistAY <= tol ) |> sapply(isTRUE)

  gb.bak$tra <- Rle(gb$tra)
  gb.bak$trs <- Rle(gb$trs)
  gb.bak
}

#' Show translocations and their flanking blocks.
#'
#' @param gb A [`GBreaks`] object.
#' @param rename Replace range names by their numeric order before subsetting.
#' @param detect Run again `flagTranslocations()` if `TRUE`, else reuse the `tra`
#'        flag or fail if it is absent.
#'
#' @family Translocation functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' part of a translocation triplet have been discarded.  If the object was missing
#' the `tra` metadata column, return the object after discarding all of its
#' ranges.
#'
#' @examples
#' showTranslocations(exampleTranslocation)
#'
#' @export

showTranslocations <- function(gb, rename = TRUE, detect = TRUE) {
  if(isTRUE(detect)) gb <- flagTranslocations(gb)
  if (is.null(gb$tra)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  flagPos <- which(gb$tra) + 1
  flagContext <- c(flagPos -1 , flagPos, flagPos + 1) |> unique() |> sort()
  gb[flagContext]
}

#' Extract central blocks in translocations
#'
#' @param gb A [`GBreaks`] object.
#' @param rename Replace range names by their numeric order before subsetting.
#' @param remove Filter out instead of filtering in.
#' @param detect Run again `flagTranslocations()` if `TRUE`, else reuse the `tra`
#'        flag or fail if it is absent.
#'
#' @family Translocation functions
#'
#' @returns Returns the `GBreaks` object in which all ranges that are not
#' the central part of an inversion triplet have been discarded.  If the object
#' was missing the `tra` metadata column, return the object after discarding all
#' of its ranges.
#'
#' @examples
#' filterTranslocations(exampleTranslocation)
#' filterTranslocations(remove = TRUE, exampleTranslocation)
#'
#' @export

filterTranslocations <- function(gb, rename = TRUE, remove = FALSE, detect = TRUE) {
  if(isTRUE(detect)) gb <- flagTranslocations(gb)
  if (is.null(gb$tra)) return(gb[0])
  if (isTRUE(rename))
    names(gb) <- seq_along(gb)
  flagPos <- which(gb$tra) + 1
  if (isTRUE(remove)) {
    # gb[-integer(0)] would return an empty object.
    if (length(flagPos) == 0) return(gb)
    flagPos <- -(flagPos)
  }
  gb[flagPos]
}

#' Remove translocations
#'
#' @param gb A [`GBreaks`] object.
#' @param detect Run again `flagTranslocations()` if `TRUE`, else reuse the `tra`
#'        flag or fail if it is absent.
#'
#' @family Translocation functions
#' @family modifier functions
#'
#' @returns Returns the `GBreaks` object after removing the ranges that are
#' translocated and coalescing the resulting collinear ranges.
#'
#' @examples
#' removeTranslocations(exampleTranslocation)
#'
#' @export

removeTranslocations <- function(gb, detect = TRUE) {
  if (length(gb) < 3) return(gb)
  if (isTRUE(detect)) gb <- flagTranslocations(gb)
  if (is.null(gb$tra)) return(gb[0])
  traPos <- which(gb$tra) + 1
  if (length(traPos) == 0) return(gb)
  coalesce_contigs(gb[-traPos])
}
