#' Sequence name with strand information
#'
#' Extract sequence names and pastes strand information to it.
#'
#' @param gb A `GBreaks` object or a [`GRangesList`] of `GBreaks`
#'           objects.
#' @param flip Flip the strand names.
#' @param query Work on the query genome.
#'
#' @return Returns a character vector, or a list of character vectors if the
#' input was a [`GenomicRanges::GRangesList`].
#'
#' @author Charles Plessy
#'
#' @family scaffolding functions
#'
#' @examples
#' strandNames(exampleColinear)
#' strandNames(exampleColinear, query = TRUE)
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
#' @family scaffolding functions
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
