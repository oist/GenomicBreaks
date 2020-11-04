#' Coverage over Breakpoints
#'
#' The function returns a GRanges object containing the breakpoints, and metadata colummns of the associated point and averaged coverage over the breakpoint. The averaged coverage is conducted over a window, and the size of which is determined by the user. The coverage is with regards to the reference genome. Make the query genome the subject of the GRanges object to extract this information for the query genome.
#'
#' @param gr_ob GRanges object containing the pairwise alignment (or alignments of the desired genome)
#' @param cov_gr GRanges object of the point coverage of the reference (or desired) genome. Point coverage information should be in the metadata column "name"
#' @param win an even number indicating the range overage which the average coverage should be calculated.
#' @return GRanges object of breakpoints and associated point and averaged coverage information
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom S4Vectors queryHits
#' @importFrom stats na.omit

# if window is greater than possible rang (i.e extends outside of a scaffold) then just give greatest possible range, even if asymmetric

bp_coverage <- function(gr_ob, cov_gr, win){

  # Define breakpoints
  gr_starts <- start(ranges(gr_ob))
  gr_ends <- end(ranges(gr_ob))
  bps_s_e <- c(gr_starts, gr_ends)
  bps_ir <- IRanges(start = bps_s_e, end = bps_s_e)
  bps_vec_seq <- as.vector(seqnames(gr_ob))
  bps_seqs <- c(bps_vec_seq, bps_vec_seq)
  bps <- GRanges(seqnames = bps_seqs, ranges = bps_ir, seqinfo = seqinfo(gr_ob))


  # find positions of breakpoints in cov_gr
  #cov_bps <- subsetByOverlaps(cov_gr, bps)
  #where_cov_bps <- which(overlapsAny(cov_gr, cov_bps) == TRUE)
  where_bps <- findOverlaps(cov_gr, bps)
  where_cov_bps <- queryHits(where_bps)
  #where_cov_bps_int <- findOverlaps()

  # determine ranges
  starts_ini <- where_cov_bps - (win*0.5)
  ends_ini <- where_cov_bps + (win*0.5 -1)

  # construct matrix of coverages of bps (cols are coverages, rows are bps)
  #########################################################################
  # Matrix of ranges
  idx_ranges <- rep.row(starts_ini, win) + rep.col(0:(win - 1), length(starts_ini))
  covs_mat <- matrix(sapply(
    idx_ranges, function(x,y) as.numeric(y[x]), y = as.vector(cov_gr$name)), nrow = dim(idx_ranges)[1], ncol = dim(idx_ranges)[2])


  # Matrix of scaffolds
  scafs_mat <- matrix(sapply(
    idx_ranges, function(x,y) y[x], y = as.vector(seqnames(cov_gr))), nrow = dim(idx_ranges)[1], ncol = dim(idx_ranges)[2])

  # Correct ranges matrix
  bp_scafs <- scafs_mat[(win*0.5),]
  bp_point_cov <- as.numeric(covs_mat[(win*0.5),])
  bp_scafs_rep <- rep.row(bp_scafs, dim(idx_ranges)[1])

  # Construct final coverage matrix
  match_scafs <- scafs_mat == bp_scafs_rep
  covs_mat_fin <- covs_mat
  covs_mat_fin[!match_scafs] <- NA

  #########################################################################

  # vector of averages
  mean_cov <- colMeans(na.omit(covs_mat_fin)) #apply(covs_mat_fin, 2, function(x) mean(na.omit(x)))

  # add meta data columns of point and averaged coverage to bps
  bps$coverage <- mean_cov
  bps$point_coverage <- bp_point_cov

  return(bps)

}

rep.row <- function(x,n){
  matrix(rep(x,each=n), nrow=n)
}

rep.col <- function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}
