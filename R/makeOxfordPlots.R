#' Make Oxford Plots
#'
#' Takes a [`GBreaks`] object and prepares an Oxford (macrosynteny) of the
#' coordinates of the _query_ ranges against the _target_ ranges after
#' concatenating them.
#'
#' @param gb A `GBreaks` object
#' @param sp1Name Name of the first species (default: sp1)
#' @param sp2Name Name of the second species (default: sp2)
#' @param sp1ChrArms A `GBreaks` object of chromosome arms in sp1 genome
#' @param sp2ChrArms A `GBreaks` object of chromosome arms in sp2 genome
#' @param type The type of the plot (`point`, `line` or `none`)
#' @param size The size of the plotted dots or segments.
#' @param diag Diagonalise the plot by reordering _query_ sequence levels.
#' @param col Colour of the plotted dots or lines by `seqnames`, `strand` or
#'        `score`.
#'
#' @return Returns a `ggplot2` object that can be further modified using the
#' `+` operator.  Use `type = 'none'` to receive an object without _geom_
#' layer.
#'
#' @author Aleksandra Bliznina
#'
#' @family plot functions
#'
#' @examples
#' makeOxfordPlots(exampleTranslocation)
#' makeOxfordPlots(exampleTranslocation, type = "p")
#'
#' @import ggplot2
#' @export

makeOxfordPlots <- function (gb, sp1Name = "target", sp2Name = "query",
                             sp1ChrArms = NULL, sp2ChrArms = NULL,
                             type = c("line", "point", "none"), size = 1,
                             diag = TRUE, col = c("seqnames", "strand", "score")) {

  col  <- match.arg(col)
  type <- match.arg(type)

  if (sp1Name == "target") {
    possibleName <- unique(genome(gb))
    if (length(possibleName) == 1 & ! is.na(possibleName))
    sp1Name <- possibleName
  }

  if (sp2Name == "query") {
    possibleName <- unique(genome(gb$query))
    if (length(possibleName) == 1 & ! is.na(possibleName))
      sp2Name <- possibleName
  }

  # Discard seqlevels not in use
  seqlevels(gb)       <- seqlevelsInUse(gb)
  seqlevels(gb$query) <- seqlevelsInUse(gb$query)

  # Ensure we have sequence lengths
  gb <- forceSeqLengths(gb)

  # Merge seq levels if needed.
  mergeSeqLevelsIfMany <- function(gr, seqs, name) {
    if (length(seqlevelsInUse(gr)) == 1) return(gr)
    mergeSeqLevels(gr, seqs, name)
  }

  targetMerged <- mergeSeqLevelsIfMany(gb,       seqlevelsInUse(gb),       sp1Name)

  if(isTRUE(diag)) {
    newOrder <- orderQuerySeqLevels(GBreaks(target = targetMerged, query = gb$query))
    seqlevels(gb$query) <- seqlevels(gb$query)[newOrder]
  }

  queryMerged  <- mergeSeqLevelsIfMany(gb$query, seqlevelsInUse(gb$query), sp2Name)

  breaks <- list()
  ticks  <- list()
  labels <- list()

  ## This is needed to silence "no visible binding for global variable" NOTEs
  start <- query.start <- end <- query.end <- seqnames <- NULL
  ## This is needed to build the data frame below regardless a score exists.
  if (is.null(score(gb))) score(gb) <- NA
  ## plot main data
  p <- ggplot(data.frame(      start = start(targetMerged),
                                 end =   end(targetMerged),
                              strand = strand(targetMerged),
                         query.start = ifelse(strand (gb) == "+", start(queryMerged),   end(queryMerged)),
                           query.end = ifelse(strand (gb) == "+",   end(queryMerged), start(queryMerged)),
                            seqnames = seqnames(gb),
                               score = score(gb))) +
    aes(x = start, y = query.start, xend = end, yend = query.end)

  if (type == "point")
    p <- p + geom_point(aes_string(colour = col), shape = 20, size = size)

  if (type == "line")
    p <- p + geom_segment(aes_string(colour = col),
                          lineend = "round",
                          size    = size)

  ## add title of the plot
  p <- p + ggtitle(paste(sp1Name, "vs", sp2Name, "macrosynteny plot")) +
    theme(plot.title = element_text(hjust = 0.5))

  ## add x and y labels
  p <- p + labs(x = seqlevelsInUse(targetMerged),
                y = seqlevelsInUse(queryMerged))

  ## add breaks
  calcBreakPosition <- function(gr) {
    # list of seqlengths, starting with zero
    lengths <- c(0,seqlengths(gr)[seqlevelsInUse(gr)])
    # Return the position of the boundary of each seqlengths
    # in the context of the merged object
    lengths |> unname()|> cumsum()
  }
  calcLabelPosition <- function(breakPos) {
    # Return the position of the midpoint of each seqlengths
    # in the context of the merged object
    breakPos |> zipWithNext() |> as.data.frame() |>
      rowMeans() |> head(-1)
  }

  if (is.null(sp1ChrArms)) {
    breaks$sp1 <- calcBreakPosition(gb)
    ticks$sp1  <- calcLabelPosition(breaks$sp1)
    labels$sp1 <- seqlevelsInUse(gb)
  } else {
    breaks$sp1 <- end(sp1ChrArms)
    labels$sp1 <- as.character(round(breaks$sp1 / 10**6, 1))
  }
  p <- p +
    scale_x_continuous(expand = c(0, 0), minor_breaks = NULL,
                       breaks = breaks$sp1, labels = NULL, position = 'top',
                       sec.axis = dup_axis(breaks=ticks$sp1, labels=labels$sp1))

  if (is.null(sp2ChrArms)) {
    breaks$sp2 <- calcBreakPosition(gb$query)
    ticks$sp2  <- calcLabelPosition(breaks$sp2)
    labels$sp2 <- seqlevelsInUse(gb$query)
  }
  else {
    breaks$sp2 <- end(sp2ChrArms)
    labels$sp2 <- as.character(round(breaks$sp2 / 10**6, 1))
  }
  p <- p +
    scale_y_continuous(expand = c(0, 0), minor_breaks = NULL,
                       breaks = breaks$sp2, labels = NULL, position = 'right',
                       sec.axis = dup_axis(breaks=ticks$sp2, labels=labels$sp2)) +
    coord_fixed()

  if (col == 'seqnames')
    p <- p + theme(legend.position="none")

  p
}
