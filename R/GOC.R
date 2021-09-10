#' Gene Order Conservation
#'
#' The Gene Order Conservation (GOC) number defined by Rocha (2003) is:
#' \dQuote{\emph{the average number of orthologues for which the consecutive
#' orthologue co-occurs close by in the other genome. It varies between 0
#' (no co-occurrence) and 1 (complete gene order conservation)}}.
#'
#' @note Note that calculating GOC on whole-genome alignments is not expected
#' to produce meaningful results.  This function is more useful when comparing
#' the position of orthologues, represented in a `GBreaks` object.
#'
#' @param gb A [`GBreaks`] object.
#'
#' @param vicinity How far to search for the neighbour orthologue.
#'
#' @param debug See below.
#'
#' @returns Returns a numeric value between 0 and 1.  If `debug = TRUE`,
#' returns a copy of the `gb` object with additional columns showing details of
#' the computation.
#'
#' @family Colinearity functions
#' @family Similarity indexes
#'
#' @references Rocha, Eduardo P C. \dQuote{DNA repeats lead to the accelerated
#' loss of gene order in bacteria.} _Trends in genetics : TIG_ vol. 19,11
#' (2003): 600-3. [doi:10.1016/j.tig.2003.09.011](https://doi.org/10.1016/j.tig.2003.09.011)
#'
#' @examples
#' gb1       <- GRanges(c(A="Ref:100-200:+", B="Ref:400-500:+"))
#' gb1$query <- GRanges(c(A="Que:100-200",   B="Que:400-500"))
#' gb1
#' GOC(gb1)
#'
#' gb4       <- GRanges(c("Ref:100-200:+",   "Ref:400-500:+",   "Ref:600-700:+"))
#' gb4$query <- GRanges(c("Que:1100-1200:+", "Que:1700-1800:+", "Que:1500-1600:+"))
#' gb4
#' GOC(gb4)
#' GOC(gb4, vicinity = 1)
#'
#' inv       <- GRanges(c("XSR:101-180:+", "XSR:201-300:-",  "XSR:320-400:+"))
#' inv$query <- GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400"))
#' inv
#' GOC(inv, vicinity = 1)
#'
#' @export

GOC <- function(gb, vicinity = 4, debug = FALSE) {
  if(length(gb) == 0) stop("Can not compute GOC on zero-length objects.")

  # Re-sort the object just to be sure
  gb <- sort(gb, ignore.strand = TRUE)

  # Get relative positions of orthologues
  gb <- flagColinearAlignments(gb, details = TRUE)

  # TRUE if neighbor orthologue within vicinity
  gb$GOCi <- abs(gb$qfoll) <= vicinity

  # Debug statement to help sanity checks
  if (debug) return(gb)

  # Remove terminal blocks on the target genome as they have no following
  # neighbors.
  gb[is.na(gb$tfoll)] <- NULL

  # Transform remaining NAs into FALSE before summing.
  sum(gb$GOCi %in% TRUE) / length(gb)
}
