#' Index representing how 'syntenic' two genomes are
#'
#' A sequence feature on the _target_ genome will be considered \sQuote{syntenic}
#' if most of its alignments map to the same feature of the _query_ genome, and
#' the two genomes are considered syntenic if most of their features are
#' syntenic.
#'
#' For a given feature of the _target_ genome the number of bases aligned on
#' each of the features of the _query_ genome are computed.  The largest number
#' is taken.  To make features comparable, this number is then normalised by
#' total number of aligned bases.  Then, the computed numbers for all the
#' _target_'s features are averaged with a weighted mean, where the weight is
#' the total number of aligned bases of that feature.  Thus, a number close to
#' 1 is expected for closely related genomes.
#'
#' @note Here, the definition of \sQuote{synteny} is _sequences on the same
#' DNA strands_, or in technical terms: sequences that are on the same feature,
#' regardless of strand orientation.
#'
#' @family Similarity indexes
#'
#' @returns Returns a numeric value between 0 and 1.
#'
#' @examples
#' gb       <- GenomicRanges::GRanges(c("Ref:100-200:+",   "Ref:400-500:+",    "Ref:600-700:+"))
#' gb$query <- GenomicRanges::GRanges(c("Que:1100-1200:+", "Que2:1700-1800:+", "Que:1500-1600:+"))
#' synteny_index(gb)
#'
#' @export

synteny_index <- function(gb) {
  gbl <- split(gb, droplevels(seqnames(gb)))
  # Calculate an index for each sequence feature
  synIdx <- sapply(gbl, \(x) {
    # Largest sum of all widths
    alnOnMainHit <- split(x, seqnames(x$query)) |> width() |> sum() |> sort() |> tail(1)
    # Normalise by the sum of all widths
    alnOnMainHit / sum(width(x))
  })
  # Average by the sum of all widths
  weighted.mean(synIdx, sum(width(gbl)))
}
