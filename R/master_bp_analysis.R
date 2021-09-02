#' Master Breakpoint Analysis
#'
#' This function brings together the different processes and analyses of the GenomicBreaks package; coalescents, breakpoint coverage and breakpoint proximity to tandem repeats. Input fields can be left blank if the information is not available, and corresponding analyses will not be outputed.
#'
#' @param pair_gr GRanges object containing pairwise alignment
#' @param co_tol tolerance over which the coalescing algorithm in \code{coalesce_contigs} should bridge. If left blank, coalescing will not be performed
#' @param ref_tan tandem repeat coverage for the reference genome in GRanges format
#' @param ref_tan_tol proximity to which a breakpoint should be considered to be "near" a tandem repeat
#' @param q_tan tandem repeat coverage for the query genome in GRanges format
#' @param q_tan_tol proximity to which a breakpoint should be considered to be "near" a tandem repeat
#' @param ref_cov reference genome coverage (per base)
#' @param ref_cov_tol range overwhich to give average coverage over a breakpoint (reference genome)
#' @param q_cov query genome coverage (per base)
#' @param q_cov_tol range overwhich to give average coverage over a breakpoint (query genome)
#' @return a list containing two GRanges objects; one for the reference genome and one for the query genome. Metadata columns contain information per breakpoint, indicated by left (start) breakpoint and right (end) breakpoint
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom stats na.omit

# need to make sure that the seqnames/scaffold names are consistent across ref genome, ref_tan, ref_cov
# and similarly between query genome, q_tan, q_cov

master_bp_analysis <- function(pair_gr, co_tol, ref_tan, ref_tan_tol, q_tan, q_tan_tol, ref_cov, ref_cov_tol, q_cov, q_cov_tol){

  #### Coalescing algorithm
  # run algorithm (only if coalescing tolerance is specified, therefore default is input gr_ob)
  if (!missing(co_tol)){
  red_gr <- coalesce_contigs(pair_gr, tol = co_tol)
  } else {
    red_gr <- pair_gr
  }

  # split into ref and query
  gr_ref <- red_gr
  gr_q <- red_gr$name
  gr_ref$name <- NULL

  #### coverage info
  ## for gr_ref
  if (!missing(ref_cov) && !missing(ref_cov_tol)){
    gr_ref_cov <- bp_coverage(gr_ob = gr_ref, cov_gr = ref_cov, win = ref_cov_tol)
    gr_ref <- give_coverage(gr_object = gr_ref, cov_gr_ob = gr_ref_cov)
  }

  ## for gr_q
  if (!missing(q_cov) && !missing(q_cov_tol)){
    gr_q_cov <- bp_coverage(gr_ob = gr_q, cov_gr = q_cov, win = q_cov_tol)
    gr_q <- give_coverage(gr_object = gr_q, cov_gr_ob = gr_q_cov)
  }

  #### tandem repeat info
  ## for gr_ref
  if (!missing(ref_tan) && !missing(ref_tan_tol)){
  ref_tan_out <- tan_bp(gr_ob = gr_ref, tan = ref_tan, tol = ref_tan_tol, query_tf = FALSE)

  # match positions for metadata columns
  ref_tan_prox <- give_tandem(tan_prox = ref_tan_out, gr_object = gr_ref)
  gr_ref$L_bp_tan_prox <- ref_tan_prox[[1]]
  gr_ref$R_bp_tan_prox <- ref_tan_prox[[2]]
  }

  ## for gr_q
  if (!missing(q_tan) && !missing(q_tan_tol)){
  q_tan_out <- tan_bp(gr_ob = gr_q, tan = q_tan, tol = q_tan_tol, query_tf = FALSE)

  # match positions for metadata columns
  q_tan_prox <- give_tandem(tan_prox = q_tan_out, gr_object = gr_q)
  gr_q$L_bp_tan_prox <- q_tan_prox[[1]]
  gr_q$R_bp_tan_prox <- q_tan_prox[[2]]
  }

  #### return granges objects

  return(list(gr_ref, gr_q))
}

######### Functions ############################################

## Proximity to tandem repeats
#tan_prox is output from tan_bp
#gr_object is either reference or query granges object
give_tandem <- function(tan_prox, gr_object){
  gr_tan_starts <- vector(mode = "list", length = length(gr_object))
  gr_tan_ends <- vector(mode = "list", length = length(gr_object))
  num_scafs <- length(levels(seqnames(gr_object)))
  k <- 0 #counter
  for (i in 1:num_scafs){
    scaf_now <- levels(seqnames(gr_object))[i]
    gr_now <- gr_object[seqnames(gr_object) == scaf_now]
    how_many_scaf_now <- length(gr_now)

    tan_in <- tan_prox[[1]][seqnames(tan_prox[[1]]) == scaf_now]
    in_match_start <- match(start(ranges(tan_in)), start(ranges(gr_now)))
    in_match_end <- match(start(ranges(tan_in)), end(ranges(gr_now)))
    gr_tan_starts[na.omit(in_match_start) + k] <- "within"
    gr_tan_ends[na.omit(in_match_end) + k] <- "within"

    tan_near <- tan_prox[[2]][seqnames(tan_prox[[2]]) == scaf_now]
    near_match_start <- match(start(ranges(tan_near)), start(ranges(gr_now)))
    near_match_end <- match(start(ranges(tan_near)), end(ranges(gr_now)))
    gr_tan_starts[na.omit(near_match_start) + k] <- "near"
    gr_tan_ends[na.omit(near_match_end) + k] <- "near"

    tan_out <- tan_prox[[3]][seqnames(tan_prox[[3]]) == scaf_now]
    far_match_start <- match(start(ranges(tan_out)), start(ranges(gr_now)))
    far_match_end <- match(start(ranges(tan_out)), end(ranges(gr_now)))
    gr_tan_starts[na.omit(far_match_start) + k] <- "far"
    gr_tan_ends[na.omit(far_match_end) + k] <- "far"

    k <- k+how_many_scaf_now
  }
  return(list(gr_tan_starts, gr_tan_ends))
}


## Coverage

give_coverage <- function(gr_object, cov_gr_ob){
  gr_object$left_cov_pb <- cov_gr_ob$point_coverage[1:(length(cov_gr_ob)*0.5)]
  gr_object$left_cov_av <- cov_gr_ob$coverage[1:(length(cov_gr_ob)*0.5)]
  gr_object$right_cov_pb <- cov_gr_ob$point_coverage[(length(cov_gr_ob)*0.5 + 1) : (length(cov_gr_ob))]
  gr_object$right_cov_av <- cov_gr_ob$coverage[(length(cov_gr_ob)*0.5 + 1) : (length(cov_gr_ob))]

  return(gr_object)
}


