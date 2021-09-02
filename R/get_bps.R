#' Breakpoints
#'
#' Given a GRanges object, the function produces a GRanges object detailing the breakpoints only. The only converted data is that of the main GRanges subject, and no metadata is processed or carried through
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param direction Will return the breakpoints on `both`, `left` or `right` side(s) of the range.
#' @param stranded If `TRUE`, will assign a `+` strand to the left-side breakpoints and a `-` strand to the right-side ones.
#' @param sorted Sorts the object before returning it.
#' @return GRanges object of the breakpoints
#' @examples
#' gr <- GRanges(c("chr2", "chr2", "chr1", "chr3"), IRanges::IRanges(1:4, width=4:1))
#' get_bps(gr)
#' get_bps(gr, direction = "left")
#' get_bps(gr, stranded = TRUE)
#' get_bps(gr, direction = "right", stranded = TRUE)
#' @importFrom GenomicRanges start end GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames
#' @export

get_bps <- function(gr_ob, direction = c("both", "left", "right"), stranded = FALSE, sorted = TRUE) {
  direction <- match.arg(direction) # stops if `direction` is not `both`, `left` or `right`
  gr_starts <- flank(gr_ob, -1, start = TRUE,  ignore.strand = TRUE) # start bps
  gr_ends   <- flank(gr_ob, -1, start = FALSE, ignore.strand = TRUE) # end bps
  if (stranded) {
    strand(gr_starts) <- "+"
    strand(gr_ends)   <- "-"
  } else {
    strand(gr_starts) <- "*"
    strand(gr_ends)   <- "*"
  }
  if (direction == "both") {
    gr <- c(gr_starts, gr_ends) # concatenate start and end bps
  } else if (direction == "left") {
    gr <- gr_starts
  } else if (direction == "right") {
    gr <- gr_ends
  }
  if (sorted) gr <- sort(gr, ignore.strand = TRUE)
  granges(gr)
}
