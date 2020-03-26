#' Breakpoints
#'
#' Given a GRanges object, the function produces a GRanges object detailing the breakpoints only. The only converted data is that of the main GRanges subject, and no metadata is processed or carried through
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @return GRanges object of the breakpoints
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb

get_bps <- function(gr_ob){

  gr_starts <- start(ranges(gr_ob)) # start bps
  gr_ends <- end(ranges(gr_ob)) # end bps
  bps_s_e <- c(gr_starts, gr_ends) # concatenate start and end bps
  bps_ir <- IRanges(start = bps_s_e, end = bps_s_e) # IRanges object needed for constructing GRanges object
  bps_vec_seq <- as.vector(seqnames(gr_ob))
  bps_seqs <- c(bps_vec_seq, bps_vec_seq)
  bps <- GRanges(seqnames = bps_seqs, ranges = bps_ir, seqinfo = seqinfo(gr_ob))
}
