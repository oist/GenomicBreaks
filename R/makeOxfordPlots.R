#' Make Oxford Plots
#'
#' TBD
#'
#' @param gb A GenomicBreaks object
#' @param selChroms A list of chromosomes in target genome to plot
#' @param sp1Name Name of the first species (default: sp1)
#' @param sp2Name Name of the second species (default: sp2)
#' @param sp1ChrArms A GenomicBreaks object of chromosome arms in sp1 genome
#' @param sp2ChrArms A Genomic Breaks object of chromosome arms in sp2 genome
#' @param type The type of the plot (`point` or `line`)
#'
#' @return Returns an Oxford (macrosynteny) plot generated with `ggplot2`
#'
#' @author Aleksandra Bliznina
#'
#' @family plot functions
#'
#' @examples
#' makeOxfordPlots(exampleTranslocation)
#' makeOxfordPlots(exampleTranslocation, type = "l")
#'
#' @import ggplot2
#' @export

makeOxfordPlots <- function (gb, selChroms = NULL,
                             sp1Name = "sp1", sp2Name = "sp2",
                             sp1ChrArms = NULL, sp2ChrArms = NULL,
                             type = c("point", "line")) {


  type <- match.arg(type)

  # filter chromosomes
  if(! is.null(selChroms)){
    gb <- gb[seqnames(gb) %in% selChroms]
  }

  targetMerged <- mergeSeqLevels(gb,       seqlevelsInUse(gb), "AllMerged")
  queryMerged  <- mergeSeqLevels(gb$query, seqlevelsInUse(gb$query), "AllMerged")

  breaks <- list()
  labels <- list()

  ## plot main data
  p <- ggplot(data.frame(      start = start(targetMerged),
                                 end =   end(targetMerged),
                         query.start = ifelse(strand (gb) == "+", start(queryMerged),   end(queryMerged)),
                           query.end = ifelse(strand (gb) == "+",   end(queryMerged), start(queryMerged)),
                            seqnames = seqnames(gb))) +
    aes(x = start, y = query.start, xend = end, yend = query.end)

  if (type == "point")
    p <- p + geom_point(aes(colour = seqnames), shape = 20, size = 0.01)

  if (type == "line")
    p <- p + geom_segment(aes(colour = seqnames))

  ## add title of the plot
  p <- p + ggtitle(paste(sp1Name, "vs", sp2Name, "macrosynteny plot")) +
    theme(plot.title = element_text(hjust = 0.5))

  ## add x and y labels
  p <- p + labs(x = sp1Name, y = sp2Name)

  ## add breaks
  if (is.null(sp1ChrArms)) {
    breaks$sp1 <- unname(c(0, cumsum(head(seqlengths(gb)[seqlevelsInUse(gb)], -1))))
    labels$sp1 <- seqlevelsInUse(gb)
  } else {
    breaks$sp1 <- end(sp1ChrArms)
    labels$sp1 <- as.character(round(breaks$sp1 / 10**6, 1))
  }
  p <- p +
    scale_x_continuous(expand = c(0, 0), breaks = breaks$sp1,
                       minor_breaks = NULL, labels = labels$sp1)

  if (is.null(sp2ChrArms)) {
    if (all(is.na(seqlengths(gb$query))))
      seqlengths(gb$query) <- tapply(end(gb$query), seqnames(gb$query), max) |> as.vector()
    breaks$sp2 <- unname(c(0, cumsum(head(seqlengths(gb$query)[seqlevelsInUse(gb$query)], -1))))
    labels$sp2 <- seqlevelsInUse(gb$query)
  }
  else {
    breaks$sp2 <- end(sp2ChrArms)
    labels$sp2 <- as.character(round(breaks$sp2 / 10**6, 1))
  }
  p <- p +
    scale_y_continuous(expand = c(0, 0), breaks = breaks$sp2,
                       minor_breaks = NULL, labels = labels$sp2)

  p
}
