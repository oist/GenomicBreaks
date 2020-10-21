#' Breakpoint Associated Heatmaps for Pairwise Aligned Genomes
#'
#' This function generates a Heatmap of the specified pattern, over breakpoints aligned at the centre of the plot. The result can be piped into \code{smoothHeatmap} and then \code{plotHeatmapList} or \code{plotHeatmapMeta}
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param gen_seq This should be either a BSgenome object or a DNAStringSet object, such that the genome sequence is contained in this input
#' @param basep_range range over which heatmap characteristic is plotted. Breakpoints will be aligned at the center of this.
#' @param pat character string of desired pattern/characteristic to be plotted on heatmap
#' @param ... Pass other arguments to \code{get_bps}.
#' @return Heatmap of pattern around centred breakpoints
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps PatternHeatmap
#' @importFrom Biostrings getSeq

bp_heatmap <- function(gr_ob, gen_seq, basep_range, pat, ...){

  # Suppress warnings about overflow
  suppressWarnings(
  gr_bps <- get_bps(gr_ob, ...) + basep_range/2
  )
  # Remove the out-of-bound ranges.
  gr_bps <- gr_bps[gr_bps == trim(gr_bps)]

  # Heatmap
  hm_prep <- getSeq(gen_seq, gr_bps)
  PatternHeatmap(hm_prep, pattern = pat, c(-basep_range/2, basep_range/2))
}
