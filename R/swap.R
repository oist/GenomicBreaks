#' Swap reference and query genomes
#'
#' Produce a new object in which the information about the reference and
#' query genomes have been inversed.
#'
#' @param gr_ob Genomic ranges representing a pairwise genome alignment.
#'
#' @return Returns a new GRanges object representing the input's query ranges
#' with the input's reference ranges in the query slot.  The strand information
#' is transferred so that the query ranges are strandless and the reference
#' are stranded.
#'
#' @examples
#' gr       <- GRanges(c("chr1:100-200:+", "chr1:300-400:+", "chr1:500-600:-"))
#' gr$query <- GRanges(c("chr5:300-400",   "chr6:100-200",   "chr7:200-300"))
#'
#' @export

swap <- function(gr_ob) {
  gr               <- gr_ob$query
  gr$query         <- granges(gr_ob)
  strand(gr)       <- strand(gr_ob)
  strand(gr$query) <- "*"
  gr
}
