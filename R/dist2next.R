#' Distances to next ranges
#'
#' Calculates the distance to the next range in the reference and query genome.
#'
#' Note that because the last range has not next neighbor, the last value is set
#' to Inf arbitrarily.
#'
#' @param gr_ob A GenomicBreaks object.
#' @param ignore.strand Calculate distance for ranges on different starands.
#'
#' @return Returns the object with two extra metadata colums, `rdist` and
#' `qdist`, containing the distance to the next range in the reference and query
#' genomes respectively.
#'
#' @examples
#' gb1       <- GenomicRanges::GRanges(c("chr1:100-200:+", "chr1:400-500:+", "chr1:700-800:-"))
#' gb1$query <- GenomicRanges::GRanges(c("chr1:100-200:+", "chr1:400-500:+", "chr1:700-800:-"))
#' dist2next(gb1)
#' dist2next(gb1, ignore.strand = TRUE)
#'
#' @importFrom GenomicRanges distance
#' @importFrom utils head tail
#' @export

dist2next <- function (gr_ob, ignore.strand = FALSE) {
  d2n <- function(gr) c(distance(head(gr, -1), tail(gr, -1), ignore.strand) + 1, Inf)
  gr_ob$rdist <- d2n(gr_ob)
  gr_ob$qdist <- d2n(gr_ob$query)
  gr_ob
}
