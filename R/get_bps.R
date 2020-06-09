#' Breakpoints
#'
#' Given a GRanges object, the function produces a GRanges object detailing the breakpoints only. The only converted data is that of the main GRanges subject, and no metadata is processed or carried through
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @return GRanges object of the breakpoints
#' @export
#' gr <- GRanges(c("chr2", "chr2", "chr1", "chr3"), IRanges(1:4, width=4:1))
#' get_bps(gr)
#' @importFrom GenomicRanges start end GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqnames

get_bps <- function(gr_ob){

  gr_starts <- start(gr_ob) # start bps
  gr_ends <- end(gr_ob) # end bps
  bps_s_e <- c(gr_starts, gr_ends) # concatenate start and end bps
  bps_ir <- IRanges(start = bps_s_e, end = bps_s_e) # IRanges object needed for constructing GRanges object
  bps_vec_seq <- seqnames(gr_ob)
  bps_seqs <- c(bps_vec_seq, bps_vec_seq)
  GRanges(seqnames = bps_seqs, ranges = bps_ir, seqinfo = seqinfo(gr_ob))
}
