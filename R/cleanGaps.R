#' Gaps between ranges
#'
#' Utility function that runs `GenomicRanges::gaps()` and then cleans its output
#' by removing strandless ranges as well as the last range that represents the
#' sequence between the end of the input object and the end of its sequence
#' levels.
#'
#' @returns Returns a strandless [`GRanges`] object representing the gaps
#' between the ranges of the input.
#'
#' @param gr A [`GBreaks`] or a [`GRanges`] object.
#'
#' @note If you find replacement provided by a package that we already import,
#' please let me know in a GitHub issue or pull request.
#'
#' @family Reducing functions
#'
#' @examples
#' cleanGaps(exampleColinear)
#' GenomicRanges::gaps(exampleColinear)
#' cleanGaps(exampleInversion)
#'
#' @author Charles Plessy
#'
#' @export

cleanGaps <- function(gr) {
  # Drop all seq info, otherwise gaps() adds artificial results on unused strands
  gr.strandless <- GRanges(seqnames = seqnames(gr), ranges = ranges(gr), strand = "*")
  gps <- gaps(gr.strandless)
  gpsList <- split(gps, seqnames(gps), drop = TRUE)
  cleanList <- endoapply(gpsList, \(x) x[-1])
  cleanGaps <- unlist(cleanList)
  seqinfo(cleanGaps) <- seqinfo(gr)
  names(cleanGaps) <- NULL
  cleanGaps
}
