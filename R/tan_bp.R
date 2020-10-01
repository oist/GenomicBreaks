#' Breakpoint proximity to Tandem Repeats
#'
#' This function classifies breakpoints in terms of its proximity to tandem repeats (but can be generalised to any property whose coverage is binary).
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param tan GRanges object containing tandem repeat (or characteristic of interest) coverage
#' @param tol range in which breakpoints will be classified as "near" to a tandem repeat
#' @param query_tf logical value that, if set to TRUE, will classify bps in query coordinates as opposed to reference (default)
#' @return A list containing 3 GRanges objects of breakpoints far, near and within tandem repeats.
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb

tan_bp <- function(gr_ob, tan, tol, query_tf = FALSE){

  # extract query as GRanges object
  if(query_tf == TRUE){
    gr_q <- gr_ob$query}
  else{gr_q <- gr_ob}

  gr_bps <- get_bps(gr_q)

  # bps in tandem repeat
  bp_tan_overlap <- overlapsAny(gr_bps, tan)
  bps_in_tan <- gr_bps[bp_tan_overlap]

  # bps near tandem repeat

  start_flank <- flank(tan, tol)
  end_flank <- flank(tan, tol, start = FALSE)
  flanks <- c(start_flank, end_flank)
  flanks_red <- reduce(flanks)
  near_range <- subsetByOverlaps(flanks_red, tan, invert = TRUE)
  # checks
  #no_tan_check <- intersect(tan, near_range)
  #int1 <- intersect(near_range, flanks_red)
  #check_flank <- int1 == near_range

  bp_neartan_overlap <- overlapsAny(gr_bps, near_range)
  bps_near_tan <- gr_bps[bp_neartan_overlap]


  # all other bps

  near_in_bps <- c(bps_in_tan, bps_near_tan)
  bps_out_tan <- subsetByOverlaps(gr_bps, near_in_bps, invert = TRUE)

  list(bps_in_tan, bps_near_tan, bps_out_tan)
}
