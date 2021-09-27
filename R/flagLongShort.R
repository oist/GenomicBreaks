#' Flag long and short arms
#'
#' _Oikopleura_ genomes are special as the long and short arms of their
#' chromosomes have different properties such as `GC` or repeat content.  It
#' can be useful to know if genomic region belongs to long or a short arm.
#'
#' @param gr A `GBreaks` or a simple `GRanges` object
#'
#' @param annot A `GRanges` file containing the coordinate of arms and their
#'        nature (such as `long`, `short`, `XSR` or `YSR`) in a `Type`
#'        _metadata column_.
#'
#' @param select One method among `first`, `last` and `arbitrary`, to decide
#'        on how to handle the regions that match both arms.
#'
#' @return Returns a modified version of the object with an extra _metadata
#' column_, `Arm` in which the information from the annotation file's `Type`
#' column was transferred.  See the examples below and in the manual of
#' [`IRanges::findOverlaps`] for details. on how regions that match both arms
#' are handled.
#'
#' @family Flagging functions
#'
#' @author Charles Plessy
#'
#' @examples
#' gb       <- GRanges(c("chr1:101-180:+", "chr1:201-300:+",
#'                       "chr1:320-400:+", "chr1:501-550:+"))
#' gb$query <- GRanges(c("chrI:101-200",   "chrI:201-300",
#'                       "chrI:301-400",   "chrI:501-550"))
#' gb <- GBreaks(gb)
#' gb
#'
#' annot <- GRanges(c("chr1:1-250", "chr1:260-500"))
#' annot$Type <- c("short", "long")
#' annot
#'
#' flagLongShort(gb, annot)
#' flagLongShort(gb, annot, select = "last")
#' flagLongShort(gb, annot, select = "arbitrary")
#'
#'
#' @export

flagLongShort <- function(gr, annot, select=c("first", "last", "arbitrary")) {
  select <- match.arg(select)
  overlapHits <- findOverlaps(gr, annot, select = select)
  gr$Arm <- factor(overlapHits, labels = annot$Type[unique(overlapHits)])
  gr
}
