#' Find longest matches from target to query genome.
#'
#' Using a [`GBreaks`] object representing the alignment of a _query_ genome
#' on a _target_ genome, finds the longest match of each sequence level
#' (representing contigs, scaffolds, etc.) of the _query_ on the _target_.
#'
#' Each sequence of the _query_ is represented only once in the output, but
#' sequences of the _target_ genome can be represented multiple times if they
#' are the longest match of multiple _query_ genome sequences.  When the _target_
#' genome is more contiguous than the _query_ genome, and if there are no
#' major structural variations between them, this will reveal arrangements of
#' colinear sequences in the query genome.
#'
#' For a more compact version of the results, the output of this function can be
#' piped to `strandNames(query = TRUE)`.
#'
#' @param gb A `GBreaks` object
#' @param min.width Minimum width of a match (on the query genome) for being
#'        considered.
#' @param min.matches Discard query sequences that have fewer longest matches than
#'        `min.matches` on the target.  Default is 2, so that only results relevant
#'        to chaining genomes are kept.
#'
#' @return Returns a [`GRangesList`] object containing one [`GBreaks`] object per
#' sequence on the query genome.
#'
#' @examples
#' exampleColinear3
#' exampleColinear3 |> longestMatchesInTarget(min.width = 0, min.matches = 1)
#'
#' @author Charles Plessy
#'
#' @family scaffolding functions
#'
#' @importFrom S4Vectors endoapply
#'
#' @export

longestMatchesInTarget <- function(gb, min.width = 1e4, min.matches = 2) {
  # Score on alignment width
  score(gb) <- width(gb)
  # Remove alignments shorter than min.width
  gb <- gb[score(gb) >= min.width]
  # Split by contig/scaffold/...
  gbl <- split(gb, seqnames(gb))
  # For each contig/scaffold/... keep the alignment with the highest score
  gbl <- endoapply(gbl, function(x) x[which.max(score(x))])
  # Swap query and target and sort
  gb2 <- sort(ignore.strand=TRUE, swap(unlist(gbl)))
  # Remove unused seqlevels
  seqlevels(gb2) <- seqlevelsInUse(gb2)
  # Split by contig/scaffold/...
  gbl2 <- split(gb2, seqnames(gb2))
  # Remove query sequences with less than min.matches
  Filter(function(x) length(x) >= min.matches, gbl2)
}
