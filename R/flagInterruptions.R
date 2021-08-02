#' Flag interruptions
#'
#' Flag ranges that start a triplet where the flanking ranges are colinear if
#' the central pair is removed.
#'
#' Here is a trivial example of an interruption
#'
#' ```
#' ┌───────────────┬────────────────┬───────────────┐
#' │ XSR:101-200:+ │ Chr1:201-300:+ │ XSR:301-400:+ │  (OKI2018_I69)
#' └───────────────┴────────────────┴───────────────┘
#' ┌───────────────┬────────────────┬───────────────┐
#' │  S1:101-200:+ │  S2:201-300:+  │  S1:301-400:+ │  (OSKA2016)
#' └───────────────┴────────────────┴───────────────┘
#' ```
#'
#' @param gb A `GenomicBreaks` object.
#'
#' @return Returns the `GenomicBreaks` object with an extra `int` metadata
#' column.
#'
#' @note Note that THIS FUNCTION RETURNS `FALSE` for _inversions_ and
#' _insertions_.
#'
#' @family Flagging functions
#'
#' @examples
#' # On the plus strand
#' int       <- GenomicRanges::GRanges(c("XSR:101-200:+", "Chr1:201-300:-",  "XSR:301-400:+"))
#' int$query <- GenomicRanges::GRanges(c( "S1:101-200",   "S2:201-300",      "S1:301-400"))
#' flagInterruptions(int)
#'
#' # On the minus strand
#' int2       <- GenomicRanges::GRanges(c("XSR:301-400:-", "Chr1:201-300:+",  "XSR:101-200:-"))
#' int2$query <- GenomicRanges::GRanges(c( "S1:301-400",   "S2:201-300",      "S1:101-200"))
#' flagInterruptions(int2)
#'
#' # No interruption
#' int3       <- GenomicRanges::GRanges(c("XSR:101-200:+", "XSR:201-300:+",  "XSR:301-400:+"))
#' int3$query <- GenomicRanges::GRanges(c( "S1:101-200",   "S1:201-300",      "S1:301-400"))
#' flagInterruptions(int3)
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede
#' @export

flagInterruptions <- function (gb) {
  gb.orig <- gb # save the original object

  # The computation are made strandless because query ranges do not hold
  # strand information
  gb$tPrec <- precede(gb,       ignore.strand = TRUE) - seq_along(gb) == 2
  gb$qPrec <- precede(gb$query, ignore.strand = TRUE) - seq_along(gb) == 2
  gb$tNext <-  follow(gb,       ignore.strand = TRUE) - seq_along(gb) == 2
  gb$qNext <-  follow(gb$query, ignore.strand = TRUE) - seq_along(gb) == 2

  # Then, we check strand information
  gb$sameStrand <- strand(gb) == c(strand(tail(gb, -2)), factor(c("*", "*")))

  gb.orig$int <- ( gb$tPrec & gb$qPrec & gb$sameStrand ) |
                 ( gb$tNext & gb$qNext & gb$sameStrand )
  gb.orig$int[is.na(gb.orig$int)] <- FALSE
  gb.orig
}
