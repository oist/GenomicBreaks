#' Tandem Repeat Coverage Heatmap
#'
#' A function that produces a heatmap of the tandem coverage around breakpoints. This can be generalised to any characteristic which displays binary coverage.
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param tan_ref GRanges object containing coverage of tandem repeat in reference genome. To aquire the analysis for the query genome, reformat the GRanges object to have the query as the main subject, and replace \code{tan_ref} with the query genome tandem coverage
#' @param win window over which to observe tandem coverage over each breakpoint
#' @param lab label for the plot
#' @return The output is a coverage heatmap that can be plotted with smoothHeatmap(output), then pipe into plotHeatmapLis
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps CoverageHeatmap

tandem_coverage <- function(gr_ob, tan_ref, win, lab){

  # extract query as gr object
  #gr_q <- GRanges(gr_ob$name)

  # reference bps gr object
  ref_starts1 <- start(ranges(gr_ob)) - (win*0.5)
  ref_starts2 <- end(ranges(gr_ob)) - (win*0.5)
  ref_starts <- c(ref_starts1, ref_starts2)
  ref_ends1 <- start(ranges(gr_ob)) + (win*0.5 -1)
  ref_ends2 <- end(ranges(gr_ob)) + (win*0.5 -1)
  ref_ends <- c(ref_ends1, ref_ends2)
  ref_ir_ob <- IRanges(start = ref_starts, end = ref_ends)

  ref_seq_vec <- as.vector(seqnames(gr_ob))
  ref_seqnames <- c(ref_seq_vec, ref_seq_vec)
  ref_bps <- GRanges(seqnames = ref_seqnames, ranges = ref_ir_ob, seqinfo = seqinfo(gr_ob))

  if (FALSE){
  # query bps gr object
  q_starts1 <- start(ranges(gr_q)) - (win*0.5)
  q_starts2 <- end(ranges(gr_q)) - (win*0.5)
  q_starts <- c(q_starts1, q_starts2)
  q_ends1 <- start(ranges(gr_q)) + (win*0.5 -1)
  q_ends2 <- end(ranges(gr_q)) + (wind*0.5 -1)
  q_ends <- c(q_ends1, q_ends2)
  q_ir_ob <- IRanges(start = q_starts, end = q_ends)

  q_seq_vec <- as.vector(seqnames(gr_ob))
  q_seqnames <- c(q_seq_vec, q_seq_vec)
  q_bps <- GRanges(seqnames = q_seqnames, ranges = q_ir_ob, seqinfo = seqinfo(gr_q))
  }

  # prep tandem repeat coverage
  #seqlevels(tan_ref) <- seqlevels(gr_ob)
  #seqlengths(tan_ref) <- seqlengths(gr_ob)
  ref_trim <- ref_bps == trim(ref_bps)
  ref_bps_trim <- ref_bps[ref_trim]

  #tan_trim <- tan_ref == trim(tan_ref)
  #tan_ref_trim <- tan_ref[tan_trim]


  combined = range(c(ref_bps_trim, tan_ref))
  seqlevels(ref_bps_trim) = as.character(seqnames(combined))
  seqlevels(tan_ref) = as.character(seqnames(combined))
  seqlengths(ref_bps_trim) = end(combined)
  seqlengths(tan_ref) = end(combined)


  # heatmap
  heat_map <- CoverageHeatmap(windows = ref_bps_trim, track = tan_ref, coords = c(-win*0.5, win*0.5), label = lab)
  return(heat_map)

}
