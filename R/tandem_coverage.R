#' Tandem Repeat Coverage Heatmap
#'
#' A function that produces a heatmap of the tandem coverage around breakpoints. This can be generalised to any characteristic which displays binary coverage.
#'
#' @param gr_ob GRanges object containing pairwise alignment
#' @param tan_ref GRanges object containing coverage of tandem repeat in reference genome. To aquire the analysis for the query genome, reformat the GRanges object to have the query as the main subject, and replace \code{tan_ref} with the query genome tandem coverage
#' @param win window over which to observe tandem coverage over each breakpoint
#' @param lab label for the plot
#' @param ... Other arguments passed to \code{get_bps}`.
#' @return The output is a coverage heatmap that can be plotted with smoothHeatmap(output), then pipe into plotHeatmapLis
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps CoverageHeatmap

tandem_coverage <- function(gr_ob, tan_ref, win, lab, ...){

  # reference bps gr object
  suppressWarnings(
    ref_bps <- get_bps(gr_ob, ...) + win/2
  )

  # Remove windows that reach boundaries.
  ref_trim <- ref_bps == trim(ref_bps)
  ref_bps_trim <- ref_bps[ref_trim]

  # heatmap
  CoverageHeatmap(windows = ref_bps_trim, track = tan_ref, coords = c(-win*0.5, win*0.5), label = lab)
}
