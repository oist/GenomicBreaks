#' Sequence name with strand information
#'
#' Extract sequence names and pastes strand information to it.
#'
#' @param gb A `GenomicBreaks` object or a `GRangesList` of `GenomicBreaks`
#'           objects.
#' @param flip Flip the strand names.
#' @param query Work on the query genome.
#'
#'
#' @return Returns a character vector, or a list of character vectors if the
#' input was a `GRangesList`.
#'
#' @examples
#' gb1       <- GRanges(c(A="Ref:100-200:+", B="Ref:400-500:+"))
#' gb1$query <- GRanges(c(A="Que:100-200",   B="Que:400-500"))
#' strandNames(gb1)
#' strandNames(gb1, query = TRUE)
#'
#' @export

strandNames <- function(gb, flip = FALSE, query = FALSE) {
  if (inherits(gb, "GRangesList")) {
    # Should better use the class system...
    return( sapply(gb, strandNames, flip, query))
  }
  if (isTRUE(query)) {
    sn <- paste0(seqnames(gb$query), strand(gb))
  } else {
    sn <- paste0(seqnames(gb), strand(gb))
  }
  if (isTRUE(flip)) {
    sn <- flipStrandNames(sn)
  }
  sn
}

#' Reverse order and flip strand
#'
#' Takes a character vector of "strand names", reverses the order of its
#' elements and flips the strand.
#'
#' @return Returns the same class as its input.
#'
#' @param x A vector of strand names.
#'
#' @export

flipStrandNames <- function (x) {
  UseMethod("flipStrandNames", x)
}

#' @export

flipStrandNames.default <- function (x) {
  x <- rev(x)
  x <- sub("-", "*", x)
  x <- sub("\\+", "-", x)
  x <- sub("\\*", "+", x)
  x
}

#' @export

flipStrandNames.list <- function (x) {
  lapply(x, flipStrandNames)
}
