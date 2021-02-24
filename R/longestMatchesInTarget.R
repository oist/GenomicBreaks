#' Find longest matches from target to query genome.
#'
#' Using a `GenomicBreaks` object representing the alignment of a query genome
#' on a target genome, finds the longest match of each sequence (representing
#' contigs, scaffolds, etc.) of the query genome on the target genome.  Thus,
#' each sequence of the query genome is represented only once in total, but
#' sequences of the target genome can be represented multiple times if they
#' are the longest match of multiple query genome sequences.  When the target
#' genome is more contiguous than the query genome, and if there are no
#' major structural variations between them, this will reveal arrangements of
#' colinear sequences in the query genome.
#'
#' @param gb A GenomicBreaks object
#' @param min.width Minimum width of a match (on the query genome) for being
#'        considered.
#' @param min.matches Discard query sequences that have fewer longest matches than
#'        `min.matches` on the target.  Default is 2, so that only results relevant
#'        to chaining genomes are kept.
#'
#' @return Returns a GRangesList object containing one GenomicBreak object per
#' sequence on the query genome.
#'
#' @details For a more compact version of the results, the output of this
#' function can be passed to the `queryStrandNames` function.
#'
#' @examples
#' \dontrun{
#' library("BSgenome.Odioica.local.OSKA2016")
#' library("BSgenome.Odioica.local.OKI2018.I69")
#' gb <- load_genomic_breaks(
#'   system.file("extdata/OSKA2016__I69-5.gff3.gz", package="GenomicBreaks"),
#'   target = OSKA2016,
#'   query = OKI2018_I69)
#' }
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
  seqlevels(gb2) <- levels(droplevels(seqnames(gb2)))
  # Split by contig/scaffold/...
  gbl2 <- split(gb2, seqnames(gb2))
  # Remove query sequences with less than min.matches
  Filter(function(x) length(x) >= min.matches, gbl2)
}
