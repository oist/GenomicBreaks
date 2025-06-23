#' Ordering permutation for _query_ relative to _target_
#'
#' Computes a permutation which rearranges the sequence levels of the _query_
#' genome so that it reflects the order in which we see the matches on a
#' pairwise plot.
#'
#' This is done by computing the average midpoint position of the _query_ ranges
#' on the _target_ genome for each seqlevel, matched by the width of the
#' _target_ ranges, so that long matches have more importance.  This avoids
#' spurious ordering due to short matches in the subtelomeric regions.
#'
#' @note The order only makes sense relative to a single sequence level of the
#' _target_ genome, so the function will stop with error if there was more
#' than one.
#'
#' @param gb A [`GBreaks`] object that has only one sequence level in use on the
#' _target_ genome.
#' @param DF A [`S4Vectors::DataFrame`] object representing the _target_ genome.
#' @param gr A [`GenomicRanges::GRanges`] object representing the _query_ genome.
#'
#' @returns Returns an integer vector of order permutations for the sequence
#' levels of the _query_ genome.
#'
#' @family seqlevel functions
#'
#' @author Charles Plessy
#'
#' @examples
#' gb       <- GRanges(c("chr1:101-180:+", "chr1:201-300:+",  "chr1:320-400:+"))
#' gb$query <- GRanges(c( "cgt8:1-100",      "ctg3:1-100",    "ctg5:1-100"))
#' gb <- GBreaks(gb)
#' # Sort alphabetically as if they were inherited from a BSgenome object
#' seqlevels(gb$query) <- c("ctg3", "ctg5", "cgt8")
#' seqlevels(gb$query)
#' # Sort by match positions on the target genome.
#' seqlevels(gb$query) <- seqlevels(gb$query)[orderQuerySeqLevels(gb)]
#' seqlevels(gb$query)
#'
#' @export

orderQuerySeqLevels <- function(gb) {
  seqlevels(gb)       <- seqlevelsInUse(gb)
  stopifnot (length(seqlevels(gb)) == 1) # Not ready for full objects
  seqlevels(gb$query) <- seqlevelsInUse(gb$query)
  gb$weight <- width(gb)
  gb$midpoint <- round(start(gb), end(gb))
  grl <- split(gb, seqnames(gb$query))
  lapply(grl, \(gr) weighted.mean(gr$midpoint, gr$weight)) |> unlist() |> order()
}

#' @rdname orderQuerySeqLevels
#' @export

orderQuerySeqLevels_DF_GR <- function(DF, gr) {
  stopifnot (length(unique(DF$seqnames)) == 1) # Not ready for full objects
  seqlevels(gr) <- seqlevelsInUse(gr)
  DF$weight <- DF$end - DF$start + 1
  DF$midpoint <- round(DF$start, DF$end)
  DFL <- split(DF, seqnames(gr))
  lapply(DFL, \(DF) weighted.mean(DF$midpoint, DF$weight)) |> unlist() |> order()
}
