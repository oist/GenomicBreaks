#' Breakpoint Pairwise Comparison using Reference Genome as Coordinate System
#'
#' This function takes in two pairwise alignments, with the same reference genome. A heatmap is produced, where the breakpoints of query1 are the centre point reference, and the breakpoints of query2 that fall within the window are projected onto it. Format of the pairwise alignments are expected to be GRanges objects, with information about the query genome in the metadata column.
#'
#' @param gr_ref_q1 pairwise alignment between the reference genome and query genome 1
#' @param gr_ref_q2 is the pw alignment between the reference genome and query genome 2
#' @param win range around query1 breakpoints of which query2 breakpoints are projected on to. Should be even number
#' @param lab string which will be the label on the graph
#' @return The output is a coverage heatmap that can be plotted with smoothHeatmap(output), then pipe into plotHeatmapList
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom schoolmath is.even
#' @importFrom heatmaps CoverageHeatmap

# We construct two GRanges objects to feed into CoverageHeatMaps

bp_pair_analysis <- function(gr_ref_q1, gr_ref_q2, win, lab){

  # Failsafes
  if(!(is.even(win))){stop("win should be an even number")}
  if(!(is.character(lab) && (length(lab) == 1))){stop("lab must be a string")}

  ## Construct g1_rco (1st GRanges object)
  # IRanges object
  q1_starts1 <- start(ranges(gr_ref_q1)) - (win*0.5)
  q1_starts2 <- end(ranges(gr_ref_q1)) - (win*0.5)
  q1_starts <- c(q1_starts1, q1_starts2)
  q1_ends1 <- start(ranges(gr_ref_q1)) + (win*0.5 -1)
  q1_ends2 <- end(ranges(gr_ref_q1)) + (win*0.5 -1)
  q1_ends <- c(q1_ends1, q1_ends2)
  q1_ir_ob <- IRanges(start = q1_starts, end = q1_ends)
  # GR object
  g1_vec_seq <- as.vector(seqnames(gr_ref_q1))
  gr1_seqnames <- c(g1_vec_seq, g1_vec_seq)
  g1_rco <- GRanges(seqnames = gr1_seqnames, ranges = q1_ir_ob, seqinfo = seqinfo(gr_ref_q1))


  ## Construct g2_rco (2nd GRanges object)
  # IRanges object
  q2_starts <- start(ranges(gr_ref_q2))
  q2_ends <- end(ranges(gr_ref_q2))
  q2_bps <- c(q2_starts, q2_ends)
  q2_ir_ob <- IRanges(start = q2_bps, end = q2_bps)
  # GR object
  g2_vec_seq <- as.vector(seqnames(gr_ref_q2))
  gr2_seqnames <- c(g2_vec_seq, g2_vec_seq)
  g2_rco <- GRanges(seqnames = gr2_seqnames, ranges = q2_ir_ob, seqinfo = seqinfo(gr_ref_q2))

  # Correcting seqlevels/lengths
  g1_i <- g1_rco == trim(g1_rco)
  g1_rco <- g1_rco[g1_i]
  g2_i <- g2_rco == trim(g2_rco)
  g2_rco <- g2_rco[g2_i]

  combined= range(c(g1_rco,g2_rco))

  seqlevels(g1_rco) = as.character(seqnames(combined))
  seqlevels(g2_rco) = as.character(seqnames(combined))
  seqlengths(g1_rco) = end(combined)
  seqlengths(g2_rco) = end(combined)

  # Heatmap
  heat_map <- CoverageHeatmap(windows = g1_rco, track = g2_rco, label = lab, coords = c(-win*0.5, win*0.5))
  return(heat_map)

}
