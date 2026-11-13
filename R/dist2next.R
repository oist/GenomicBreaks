#' Distances to next ranges
#'
#' Calculates the distance to the next range in the _target_ and _query_ genome.
#'
#' The distance is defined by the [`GenomicRanges::distance`] function: it is
#' the number of bases between the two ranges, so it is `0` for adjacent
#' ranges, and also for overlapping ranges.  It is `NA` for ranges on different
#' sequences, or on different strands unless `ignore.strand = TRUE`.  Note
#' that because the last `step` ranges (all of them if there are fewer) have
#' no next neighbor at that distance, their value is set to `Inf` arbitrarily.
#'
#' @param x A [`GenomicRanges::GRanges`] or a [`GBreaks`] object.
#' @param step Distance to the \ifelse{html}{\out{n<sup>th</sup>}}{\eqn{2n^th}}
#'        block (default: first).
#' @param ignore.strand Calculate distance for ranges on different strands.
#'
#' @return For `GRanges`, returns the object with one extra metadata colums,
#' `dist`, and for `GBreaks` two extra columns `tdist` and `qdist`, containing
#' the distance to the next range in the reference and query genomes
#' respectively.
#'
#' @family Colinearity functions
#'
#' @examples
#' dist2next(exampleInversion)
#' dist2next(granges(exampleInversion))
#' dist2next(exampleInversion, ignore.strand = TRUE)
#' dist2next(exampleInversion - 20, ignore.strand = TRUE)
#' dist2next(exampleInversion, 2)
#'
#' @importFrom GenomicRanges distance
#' @importFrom utils head tail
#' @export


setGeneric("dist2next", function(x, step = 1, ignore.strand = FALSE) standardGeneric("dist2next"))

# Number of bases between each range and the range `step` positions after it;
# Inf for the last `step` ranges (all of them if there are fewer), which have
# no such range.
.dist2next <- function(x, step, ignore.strand) {
  c( distance( head(x, -step)
             , tail(x, -step)
             , ignore.strand = ignore.strand)
   , rep(Inf, min(step, length(x)))
   )
}

dist2next_GRanges <- function(x, step = 1, ignore.strand = FALSE) {
  x$dist <- .dist2next(x, step, ignore.strand)
  x
}

dist2next_GBreaks <- function(x, step = 1, ignore.strand = FALSE) {
  # Both columns are set in one assignment because each modification of the
  # object runs its validity checks.
  mcols(x)[c("tdist", "qdist")] <-
    DataFrame( tdist = .dist2next(x,       step, ignore.strand)
             , qdist = .dist2next(x$query, step, ignore.strand))
  x
}

#' @rdname dist2next
#' @export

setMethod("dist2next", "GRanges", dist2next_GRanges)

#' @rdname dist2next
#' @export

setMethod("dist2next", "GBreaks", dist2next_GBreaks)

# Versions of GenomicRanges::precede() and follow() for ranges that may have a
# width of zero, such as mapped regions between adjacent alignments (see
# bridgeRegions()).  precede() and follow() report a zero-width range as
# preceding and following itself, and break ties by index between a zero-width
# range and a range that starts (or ends) at the same position.  Here, each
# zero-width range is placed on a point between the two bases that surround
# it, in a coordinate system where coordinates are replaced by their rank (to
# avoid integer overflows) and doubled, which preserves their order.  Objects
# without zero-width ranges are passed unchanged to precede() and follow().

zeroWidthAsPoints <- function(x) {
  zw <- width(x) == 0
  if (!any(zw)) return(x)
  s <- start(x)
  e <- end(x)
  u <- sort(unique(c(s, e)))
  rs <- match(s, u)
  re <- match(e, u)
  GRanges(seqnames(x),
          IRanges(start = ifelse(zw, 2L * rs - 1L, 2L * rs),
                  end   = ifelse(zw, 2L * rs - 1L, 2L * re)),
          strand = strand(x))
}

precedeZW <- function(x, ...) precede(zeroWidthAsPoints(x), ...)
followZW  <- function(x, ...) follow( zeroWidthAsPoints(x), ...)
