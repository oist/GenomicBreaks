#' Flag long and short arms
#'
#' _Oikopleura_ genomes are special as the long and short arms of their
#' chromosomes have different properties such as `GC` or repeat content.  It
#' can be useful to know if genomic region belongs to long or a short arm.
#'
#' @param gr A [`GBreaks`] or a simple [`GRanges`] object
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
#' annot <- GRanges(c("chrA:1-140", "chrA:150-500"), Type = c("short", "long"))
#'
#' flagLongShort(exampleColinear, annot)
#' flagLongShort(exampleColinear, annot, select = "last")
#' flagLongShort(exampleColinear, annot, select = "arbitrary")
#'
#' @export

flagLongShort <- function(gr, annot, select=c("first", "last", "arbitrary")) {
  select <- match.arg(select)
  overlapHits <- findOverlaps(gr, annot, select = select)
  annotsInUse <- overlapHits |> unique() |> Filter(f = \(x) !is.na(x))
  gr$Arm <- factor(overlapHits, labels = annot$Type[annotsInUse])
  gr
}
