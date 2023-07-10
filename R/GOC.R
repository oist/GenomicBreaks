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
#' exampleColinear
#' GOC(exampleColinear)
#'
#' exampleTranslocation
#' GOC(exampleTranslocation)
#' GOC(exampleTranslocation, v=1)
#'
#' # GOC computation is strandless
#' exampleInversion
#' GOC(exampleInversion)
#'
#' @export

GOC <- function(gb, vicinity = 4, debug = FALSE) {
  if(length(gb) == 0) return(numeric(0))
  if(length(gb) == 1) return(    1     )

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
