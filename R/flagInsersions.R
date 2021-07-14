#' Flag insersions
#'
#' Flag ranges that start a triplet where the reference range is colinear and
#' the query range is made of two blocks that would be colinear if they were
#' not separated by an block on a different molecule.
#'
#' Here is a trivial example of an insersion.
#'
#' ```
#' ┌───────────────┬───────────────┬───────────────┐
#' │ XSR:101-200:+ │ XSR:201-300:+ │ XSR:301-400:+ │  (OKI2018_I69)
#' └───────────────┴───────────────┴───────────────┘
#' ┌───────────────┬───────────────┬───────────────┐
#' │  S1:101-200:+ │  S2:801-900:+ │  S1:201-300:+ │  (OSKA2016)
#' └───────────────┴───────────────┴───────────────┘
#' ```
#'
#' @param gr_ob A `GenomicBreaks` object.
#' @param tol Tolerance window for the distance between two ranges.
#'
#' @return Returns the `GenomicBreaks` object with an extra `ins` metadata
#' column.
#'
#' @family Flagging functions
#'
#' @examples
#' ins       <- GenomicRanges::GRanges(c("XSR:101-200:+", "XSR:201-300:+",  "XSR:301-400:+"))
#' ins$query <- GenomicRanges::GRanges(c( "S1:101-200",   "S2:801-900",     "S1:201-300"   ))
#' flagInsersions(ins)
#'
#' @include dist2next.R
#' @importFrom GenomicRanges precede
#' @export

flagInsersions <- function (gr_ob, tol = Inf) {
  gr <- gr_ob # save the original object
  gr_ob <- dist2next(gr_ob)
  gr_ob$tdist2 <- c(gr_ob$tdist[-1], Inf)
  gr_ob$qdist2 <- dist2next(gr_ob, step = 2)$qdist # Important difference with flagIn*s*ertions.
  gr_ob$twoNext <- precede(gr_ob$query, tail(gr_ob$query, -2), ignore.strand = T) - seq_along(gr_ob) == 0
  gr$ins <- gr_ob$twoNext          &
            gr_ob$tdist  < tol + 1 &
            gr_ob$tdist2 < tol + 1 &
            gr_ob$qdist2 < tol + 1
  if (any(is.na(gr$ins)))
    gr[is.na(gr$ins)]$ins <- FALSE
  gr
}
