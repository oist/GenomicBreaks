#' Make an inversion
#'
#' @param gr [`GenomicRanges`] object
#' @param dist Distribution function for new breakpoints on `gr`'s sequence levels.
#'
#'
#' @examples
#' gr <- GRanges("chr1:1-1000000:+")
#' gr2 <- GRanges(c("chr1:1-1000000:+", "chr2:1-1000000:+"))
#' gr3 <- GRanges(c("chr1:1-1000000:+", "chr1:1000000-2000000:+"))
#' gr4 <- GRanges(c("chr1:1-1000000:+", "chr1:1000000-2000000:-"))
#' gr5 <- GRanges(c("chr1:1-1000000:-"))
#'
#' @importFrom stats runif
#' @export


makeInversion.chr <- function(gr, dist=runif) {
  if (! length(seqlevels(gr)) == 1 )
    stop ("The GRanges object must contain only one sequence level")
  startPos <- start(range(gr, ignore.strand=TRUE))
  endPos   <- end(range(gr, ignore.strand=TRUE))
  breaks <-sort(round(runif(2, min=startPos, max=endPos)))
  break1 <- c( GRanges(seqlevels(gr), IRanges(start=breaks[1]), strand = "+"),
               GRanges(seqlevels(gr), IRanges(start=breaks[1]), strand = "-"))
  break2 <- c( GRanges(seqlevels(gr), IRanges(start=breaks[2]), strand = "+"),
               GRanges(seqlevels(gr), IRanges(start=breaks[2]), strand = "-"))
  broken <- setdiff(gr, break1)
  broken <- setdiff(broken, break2)
  broken$Part                            <- "Inv"
  broken[end(broken)   < breaks[1]]$Part <- "Before"
  broken[start(broken) > breaks[2]]$Part <- "After"
  broken$Part <- factor(broken$Part, levels=c("After", "Inv", "Before"))
  grl <- split(broken, broken$Part)

  invertPos <- function(pos, start, end) end - pos + start

  grl$Inv <-
    GRanges(
      seqnames(grl$Inv),
      IRanges(
        start = invertPos(end(grl$Inv),   breaks[1], breaks[2]),
        end   = invertPos(start(grl$Inv), breaks[1], breaks[2])),
      strand = strand(invertStrand(grl$Inv))
    )
  sort(granges(unname(unlist(grl))), ignore.strand=TRUE)
}

# w100 <- gr
# for (n in 1:100) w100 <- makeInversion.chr(w100)
# library(ggplot2)
# ggplot(data.frame(width = width(w100))) +
#   aes(width) +
#   geom_histogram() +
#   scale_x_log10()
#
# w1000 <- gr
# for (n in 1:1000) w1000 <- makeInversion.chr(w1000)
# library(ggplot2)
# ggplot(data.frame(width = width(w1000))) +
#   aes(width) +
#   geom_histogram() +
#   scale_x_log10()
#
# w10000 <- gr
# for (n in 1:10000) w10000 <- makeInversion.chr(w10000)
# library(ggplot2)
# ggplot(data.frame(width = width(w10000))) +
#   aes(width) +
#   geom_histogram() +
#   scale_x_log10()

