#' Breakpoint Associated Heatmaps for Pairwise Aligned Genomes
#'
#' This function generates a Heatmap of the specified pattern, over breakpoints aligned at the centre of the plot. The result can be piped into \code{smoothHeatmap} and then \code{plotHeatmapList} or \code{plotHeatmapMeta}
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param gen_seq This should be either a BSgenome object or a DNAStringSet object, such that the genome sequence is contained in this input
#' @param basep_range range over which heatmap characteristic is plotted. Breakpoints will be aligned at the center of this.
#' @param pat character string of desired pattern/characteristic to be plotted on heatmap
#' @return Heatmap of pattern around centred breakpoints
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps PatternHeatmap
#' @importFrom Biostrings getSeq

bp_heatmap <- function(gr_ob, gen_seq, basep_range, pat){

  # Remove breakpoints that occur too close to the edges of scaffolds
  where_seqs <- match(as.vector(seqnames(gr_ob)), seqlevels(gr_ob))
  all_seq_len <- as.numeric(seqlengths(gr_ob))[where_seqs]

  ind_start1 <- start(gr_ob) < basep_range*0.5
  ind_start2 <- (start(gr_ob) + (basep_range*0.5 -1)) > all_seq_len
  ind_end1 <- (end(gr_ob) + (basep_range*0.5 -1)) > all_seq_len
  ind_end2 <- end(gr_ob) < basep_range*0.5

  mas_ind <- ind_start1 + ind_start2 + ind_end1 + ind_end2
  fin_ind <- mas_ind == 0

  red_gr <- gr_ob[fin_ind]

  # Construct new GRanges object from red_gr

  gr_bps <- get_bps(red_gr) + basep_range*0.5

  # Heatmap
  hm_prep <- getSeq(gen_seq, gr_bps)
  hm_fin <- PatternHeatmap(hm_prep, pattern = pat, c(-basep_range,basep_range))

  return(hm_fin)
}
