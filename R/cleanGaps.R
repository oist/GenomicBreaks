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
#' @param gb A [`GBreaks`] object.
#'
#' @note If you find replacement provided by a package that we already import,
#' please let me know in a GitHub issue or pull request.
#'
#' @examples
#' cleanGaps(exampleColinear)
#' GenomicRanges::gaps(exampleColinear)
#' cleanGaps(exampleInversion)
#'
#' @author Charles Plessy
#'
#' @export

cleanGaps <- function(gb) {
  # Drop all seq info, otherwise gaps() adds artificial results on unused strands
  gb <- GRanges(seqnames = seqnames(gb), ranges = ranges(gb), strand = "*")
  gps <- gaps(gb)
  gpsList <- split(gps, seqnames(gps), drop = TRUE)
  cleanList <- endoapply(gpsList, \(x) x[-1])
  unlist(cleanList)
}
