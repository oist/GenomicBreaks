#' Distances to next ranges
#'
#' Calculates the distance to the next range in the reference and query genome.
#'
#' Note that because the last range has not next neighbor, the last value is set
#' to Inf arbitrarily.
#'
#' @param gr_ob A GenomicBreaks object.
#' @param step Distance to the n^th^ block (default: first).
#' @param ignore.strand Calculate distance for ranges on different starands.
#'
#' @return Returns the object with two extra metadata colums, `tdist` and
#' `qdist`, containing the distance to the next range in the reference and query
#' genomes respectively.
#'
#' @examples
#' gb1       <- GenomicRanges::GRanges(c("chr1:100-200:+", "chr1:201-300:+", "chr1:400-500:+", "chr1:700-800:-"))
#' gb1$query <- GenomicRanges::GRanges(c("chr1:100-200:+", "chr1:201-300:+", "chr1:400-500:+", "chr1:700-800:-"))
#' dist2next(gb1)
#' dist2next(gb1, ignore.strand = TRUE)
#' dist2next(gb1, 2)
#'
#' @importFrom GenomicRanges distance
#' @importFrom utils head tail
#' @export

dist2next <- function (gr_ob, step = 1, ignore.strand = FALSE) {
  d2n <- function(gr, step) {
    c( distance( head(gr, -step)
               , tail(gr, -step)
               , ignore.strand) + 1  # Why is there +1 here ??
     , rep(Inf, step)
     )
  }
  gr_ob$tdist <- d2n(gr_ob      , step)
  gr_ob$qdist <- d2n(gr_ob$query, step)
  gr_ob
}
