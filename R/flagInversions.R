#' Flag inversions
#'
#' Flag ranges that start a triplet that would be colinear if the central pair
#' were aligned to the opposite strand.
#'
#' Here is a trivial example of an inversion.
#'
#' ```
#' ┌───────────────┬───────────────┬───────────────┐
#' │ XSR:101-200:+ │ XSR:201-300:- │ XSR:301-400:+ │  (OKI2018_I69)
#' └───────────────┴───────────────┴───────────────┘
#' ┌───────────────┬───────────────┬───────────────┐
#' │  S1:101-200:+ │  S1:201-300:+ │  S1:301-400:+ │  (OSKA2016)
#' └───────────────┴───────────────┴───────────────┘
#' ```
#'
#' @param gr_ob A `GenomicBreaks` object.
#' @param tol Tolerance window for the distance between two ranges.
#'
#' @return Returns the `GenomicBreaks` object with an extra `inv` metadata
#' column.
#'
#' @examples
#' inv       <- GenomicRanges::GRanges(c("XSR:101-200:+", "XSR:201-300:-",  "XSR:301-400:+"))
#' inv$query <- GenomicRanges::GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400"))
#' flagInversions(inv)
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede
#' @export

flagInversions <- function (gr_ob, tol = Inf) {
  gr <- gr_ob # save the original object
  gr_ob <- dist2next(gr_ob, ignore.strand = TRUE)
  gr_ob$rdist2 <- c(gr_ob$rdist[-1], Inf)
  gr_ob$qdist2 <- c(gr_ob$qdist[-1], Inf)
  gr_ob$oneNext <- precede(gr_ob$query, tail(gr_ob$query, -1), ignore.strand = T) - seq_along(gr_ob) == 0
  gr_ob$twoNext <- c(gr_ob$oneNext[-1], NA)
  gr_ob$oneDiffStrand <- strand(gr_ob) != c(strand(tail(gr_ob, -1)), factor("*"))
  gr_ob$twoSameStrand <- strand(gr_ob) == c(strand(tail(gr_ob, -2)), factor(c("*", "*")))
  gr$inv <- gr_ob$oneNext          &
            gr_ob$twoNext          &
            gr_ob$oneDiffStrand    &
            gr_ob$twoSameStrand    &
            gr_ob$rdist  < tol + 1 &   # We need `< tol + 1`
            gr_ob$rdist2 < tol + 1 &   # and not `<= tol`
            gr_ob$qdist  < tol + 1 &   # Because of comparisons to
            gr_ob$qdist2 < tol + 1     # infinite distance returned by dist2next.
  if (any(is.na(gr$inv)))
    gr[is.na(gr$inv)]$inv <- FALSE
  gr
}