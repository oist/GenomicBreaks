#' Reverse genomic ranges
#'
#' Reverse genomic ranges by flipping the strand and moving the origin of the
#' coordinate system to the opposite side.
#'
#' @returns Returns the modified `GBreaks` or `GRanges` object.
#'
#' @param x A `GBreaks` or a `GRanges` object
#' @param query On `GBreaks` objects, operate on othe _query_ or the _target_.
#' @param ... Additional arguments (ignored)
#'
#' @author Charles Plessy
#'
#' @examples
#' gb <-       GRanges(c("chr1:101-180:+", "chr1:201-300:-",  "chr2:120-400:+"))
#' gb$query <- GRanges(c("chrA:101-180",   "chrB:201-300",    "chrA:120-400"))
#' seqlengths(gb)       <- c(500, 1000)
#' seqlengths(gb$query) <- c(500, 1000)
#' gb <- GBreaks(gb)
#'
#' gb
#' reverse(gb)
#' reverse(gb, query = TRUE)
#' # Do not run reverse on the query directly as it will miss strand information
#' reverse(gb$query)
#'
#' @family modifier functions
#' @seealso See also the [`IRanges::reverse`] function.
#'
#' @importFrom IRanges reverse
#' @importMethodsFrom IRanges reverse

setGeneric("reverse", function(x, ...) standardGeneric("reverse"))

reverse_GRanges <- function(x, ...) {
  if (any(is.na(seqlengths(x))))
    stop("Can not reverse objets where `seqlengths` are not defined.")
  # Flip strand
  levels(strand(x)) <- c("-", "+", "*")

  # Opposite-side coordinates
  End <- seqnames(x)
  levels(End) <- seqlengths(x)
  End <- as.integer(as.character(End))

  # Two-step change to avoid invalid IRanges
  newStart  <- End - end(x)   + 1
  newEnd    <- End - start(x) + 1

  start(x) <- 1
  end(x)   <- End

  start(x) <- newStart
  end(x) <- newEnd

  x
}

reverse_GBreaks <- function(x, query = FALSE, ...) {
  if (isFALSE(query)) {
    gr <- reverse(granges(x))
    mcols(gr) <- mcols(x)
    x <- GBreaks(gr)
  } else {
    gr <- x$query
    strand(gr) <- strand(x)
    gr <- reverse(gr)
    strand(x) <- strand(gr)
    strand(gr) <- "*"
    x$query <- gr
  }
  x
}

#' @rdname reverse
#' @export

setMethod("reverse", "GRanges", reverse_GRanges)

#' @rdname reverse
#' @export

setMethod("reverse", "GBreaks", reverse_GBreaks)
