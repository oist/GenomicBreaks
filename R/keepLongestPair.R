#' Keep the longest pair of features
#'
#' This feature is useful when analysing genomes comprising a main chromosomes
#' and extra plasmids or smaller chromosomes.  Assuming that the main
#' chroomosome can be detected only by its length, this function extracts the
#' alignment between the pair of main chromosomes only.
#'
#' @param gb A [`GBreaks`] object
#'
#' @returns A `GBreaks` object with only one sequence feature in use per genome.
#'
#' @author Charles Plessy
#'
#' @family modifier functions
#'
#' @examples
#' exampleTranslocation
#' exampleTranslocation |> keepLongestPair()
#'
#' @export

keepLongestPair <- function(gb) {
  longestSeqFeature <- function(gr)
    seqlengths(gr) |> sort(decreasing = TRUE) |> head(1) |> names()
  gb <- forceSeqLengths(gb)
  gb[seqnames(gb)       == longestSeqFeature(gb)  &
     seqnames(gb$query) == longestSeqFeature(gb$query)]
}
