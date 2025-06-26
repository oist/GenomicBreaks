#' Breakpoint Pairwise Comparison using Reference Genome as Coordinate System
#'
#' This function takes two [`GBreaks`] objects, with the same _target_ genome.
#' A heatmap is produced, where the breakpoints of _query_ 1 are the centre
#' point reference, and the breakpoints of _query_ 2 that fall within the window
#' are projected onto it.
#'
#' @param gb1 `GBreaks` object of the alignment between the _target_ genome and _query_ genome 1.
#' @param gb2 `GBreaks` object of the alignment between the _target_ genome and _query_ genome 2.
#' @param window range around query1 breakpoints of which query2 breakpoints are projected on to. Should be even number
#' @param label string which will be the label on the graph
#'
#' @return Returns a [`heatmaps::Heatmap`] object of `pattern` that can be piped into
#' [`heatmaps::smoothHeatmap`] and then [`heatmaps::plotHeatmapList`] or
#' [`heatmaps::plotHeatmapMeta`].
#'
#' @author Charlotte West
#'
#' @family plot functions
#' @family heatmap functions
#'
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps CoverageHeatmap
#'
#' @export

bp_pair_analysis <- function(gb1, gb2, window, label){

  # Failsafes
  if( window/2 != floor(window/2) ) stop("win should be an even number")
  if(!(is.character(label) && (length(label) == 1))){stop("lab must be a string")}

  ## Construct g1_rco (1st GBreaks object)
  suppressWarnings(
  g1_rco <- get_bps(gb1) + window/2
  )
  g1_rco <- g1_rco[g1_rco == trim(g1_rco)]

  ## Construct g2_rco (2nd GBreaks object)
  suppressWarnings(
  g2_rco <- get_bps(gb2)
  )
  g2_rco <- g2_rco[g2_rco == trim(g2_rco)]

  combined <- range(c(g1_rco,g2_rco))

  seqlevels(g1_rco) = as.character(seqnames(combined))
  seqlevels(g2_rco) = as.character(seqnames(combined))
  seqlengths(g1_rco) = end(combined)
  seqlengths(g2_rco) = end(combined)

  # Heatmap
  CoverageHeatmap(windows = g1_rco, track = g2_rco, label = label, coords = c(-window*0.5, window*0.5))
}
