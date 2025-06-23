#' Get breakpoints
#'
#' Given a [`GBreaks`] or [`GenomicRanges::GRanges`] object, the function produces a `GRanges`
#' object detailing the breakpoints only.
#'
#' @param gr `GRanges` object containing pairwise alignment
#' @param direction Return the breakpoints on `both`, `left` or `right` side(s)
#'        of the range, or at the `mid`point between ranges.
#' @param stranded If `TRUE`, will assign a `+` strand to the left-side
#'        breakpoints and a `-` strand to the right-side ones.
#' @param sorted Sorts the object before returning it.
#'
#' @return `GRanges` object of the breakpoints
#'
#' @family Reducing functions
#'
#' @examples
#' get_bps(exampleInversion)
#' get_bps(exampleInversion, direction = "left")
#' get_bps(exampleInversion, stranded = TRUE)
#' get_bps(exampleInversion, direction = "right", stranded = TRUE)
#' get_bps(exampleInversion, direction = "mid")
#'
#' @importFrom GenomicRanges resize start end GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames
#'
#' @export

get_bps <- function(gr, direction = c("both", "left", "right", "mid"), stranded = FALSE, sorted = TRUE) {
  direction <- match.arg(direction)
  if (direction == "mid") {
    if(isTRUE(stranded)) stop(dQuote('direction == "mid"'), " is not compatible with ", dQuote("stranded == TRUE"), ".")
    gr <- cleanGaps(gr)
    gr <- resize(gr, 1, fix = "center")
  } else {
    gr_starts <- flank(gr, -1, start = TRUE,  ignore.strand = TRUE) # start bps
    gr_ends   <- flank(gr, -1, start = FALSE, ignore.strand = TRUE) # end bps
    if (stranded) {
      strand(gr_starts) <- "+"
      strand(gr_ends)   <- "-"
    } else {
      strand(gr_starts) <- "*"
      strand(gr_ends)   <- "*"
    }
    gr <- switch(direction,
      left  =   gr_starts,
      both  = c(gr_starts, gr_ends),
      right =              gr_ends
    )
  }
  if (sorted) gr <- sort(gr, ignore.strand = TRUE)
  granges(gr)
}
