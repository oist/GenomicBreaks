#' Distances to next ranges
#'
#' Calculates the distance to the next range in the _target_ and _query_ genome.
#'
#' The distance is defined by the [`GenomicRanges::distance`] function.  Note
#' that because the last range has not next neighbor, the last value is set
#' to `Inf` arbitrarily.
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

dist2next_GRanges <- function(x, step = 1, ignore.strand = FALSE) {
  x$dist <-
    c( distance( head(x, -step)
               , tail(x, -step)
               , ignore.strand) + 1  # Why is there +1 here ??
     , rep(Inf, step)
     )
  x
}

dist2next_GBreaks <- function(x, step = 1, ignore.strand = FALSE) {
  x$tdist <- dist2next(granges(x)      , step, ignore.strand = ignore.strand)$dist
  x$qdist <- dist2next(granges(x$query), step, ignore.strand = ignore.strand)$dist
  x$query$dist <- NULL
  x
}

#' @rdname dist2next
#' @export

setMethod("dist2next", "GRanges", dist2next_GRanges)

#' @rdname dist2next
#' @export

setMethod("dist2next", "GBreaks", dist2next_GBreaks)
