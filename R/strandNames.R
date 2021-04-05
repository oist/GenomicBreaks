#' Sequence name with strand information
#'
#' Extract sequence names and pastes strand information to it.
#'
#' @param gb A `GenomicBreaks` object of a `GRangesList` of `GenomicBreaks`
#' objects
#'
#' @return Returns a character vector, or a list of character vectors if the
#' input was a `GRangesList`.
#'
#' @examples
#' gb1       <- GenomicRanges::GRanges(c(A="Ref:100-200:+", B="Ref:400-500:+"))
#' gb1$query <- GenomicRanges::GRanges(c(A="Que:100-200",   B="Que:400-500"))
#' strandNames(gb1)
#' queryStrandNames(gb1)
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

#' @export

queryStrandNames <- function(gb, flip = FALSE) {
  strandNames(gb = gb, flip = flip, query = TRUE)
}

#' @export

flipStrandNames <- function (x, ...) {
  UseMethod("flipStrandNames", x)
}

#' Reverse order and flip strand
#'
#' Takes a character vector of "strand names", reverses the order of its
#' elements and flips the strand.
#'
#' @return Returns the same class as its input.
#'
#' @export

flipStrandNames.default <- function (sn) {
  sn <- rev(sn)
  sn <- sub("-", "*", sn)
  sn <- sub("\\+", "-", sn)
  sn <- sub("\\*", "+", sn)
  sn
}
#' @export

flipStrandNames.list <- function (l) {
  lapply(l, flipStrandNames)
}
