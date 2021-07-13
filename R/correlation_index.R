#' Index representing how correlated syntenic regions are
#'
#' A sequence feature on the _target_ genome is considered \sQuote{syntenic}
#' to the feature of the _query_ genome matched by the most of its alignments.
#' The correlation coefficient of the center of syntenic alignments is then
#' calculated.  Then, the computed numbers for all the _target_'s features are
#' averaged as a mean weighted by is the total number of aligned bases of each
#' feature.  Thus, a number close to 1 is expected for closely related genomes.
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
#' correlation_index(gb)
#'
#' @importFrom  weights wtd.cor
#'
#' @export

correlation_index <- function(gb) {
  gbl <- split(gb, droplevels(seqnames(gb)))
  # Calculate an index for each sequence feature
  corIdx <- sapply(gbl, \(x) {
    # Keep only matches to the best query
    bestQuery <- x$query |> seqnames() |> table() |> sort() |> tail(1) |> names()
    x <- x[seqnames(x$query) == bestQuery]
    if (length(x) == 1) return (1)
    # Save widths, to be used as weights later
    weights <- width(x)
    #  calculate center coordinates
    x <- resize(x, fix="center", width = 1)
    x$query <- resize(x$query, fix="center", width = 1)
    weights::wtd.cor(start(x), start(x$query), weights)[1]
  })
  # Average by the sum of all widths
  weighted.mean(corIdx, sum(width(gbl)))
}
