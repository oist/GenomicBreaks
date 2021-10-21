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
#'
#' @return Returns an Oxford (macrosynteny) plot generated with `ggplot2`
#'
#' @author Aleksandra Bliznina
#'
#' @family plot functions
#'
#' @examples
#' gb       <- GRanges(c("XSR:101-200:+", "XSR:201-300:-",  "XSR:301-400:+", "chrUn_1:401-500:+"))
#' gb$query <- GRanges(c( "S1:101-200",    "S1:201-300",     "S1:301-400", "S1:401-500"))
#' makeOxfordPlots(gb)
#'
#' @export

makeOxfordPlots <- function (gb, selChroms = NULL,
                             sp1Name = "sp1", sp2Name = "sp2",
                             sp1ChrArms = NULL, sp2ChrArms = NULL) {

  # filter chromosomes
  if(! is.null(selChroms)){
    gb <- dplyr::filter(gb, seqnames %in% selChroms)
  }

  targetMerged <- mergeSeqLevels(gb,       seqlevelsInUse(gb), "AllMerged")
  queryMerged  <- mergeSeqLevels(gb$query, seqlevelsInUse(gb$query), "AllMerged")

  breaks <- SimpleList()
  labels <- SimpleList()

  ## plot main data
  p <- ggplot(data.frame(      start = start(targetMerged),
                         query.start = start(queryMerged),
                            seqnames = seqnames(gb))) +
    aes(x = start, y = query.start) +
    geom_point(aes(colour = seqnames), shape = 20, size = 0.01)

  ## add title of the plot
  p <- p + ggtitle(paste(sp1Name, "vs", sp2Name, "macrosynteny plot")) +
    theme(plot.title = element_text(hjust = 0.5))

  ## add x and y labels
  p <- p + labs(x = sp1Name, y = sp2Name)

  ## add breaks
  if(! is.null(sp1ChrArms)){
    breaks$sp1 <- end(sp1ChrArms)
    labels$sp1 <- as.character(round(breaks$sp1 / 10**6, 1))
    p <- p +
      scale_x_continuous(expand = c(0, 0), breaks = breaks$sp1, labels = labels$sp1)
  }

  if(! is.null(sp2ChrArms)){
    breaks$sp2 <- end(sp2ChrArms)
    labels$sp2 <- as.character(round(breaks$sp2 / 10**6, 1))
    p <- p +
      scale_y_continuous(expand = c(0, 0), breaks = breaks$sp2, labels = labels$sp2)
  }
  p
}
