#' Feature coverage heatmap
#'
#' Produces a heatmap of a genomic feature around breakpoints. The feature has
#' binary coverage represented by a genomic range.
#'
#' @param gr_ob GenomicBreaks object containing pairwise alignment
#' @param feat GRanges object containing the feature in target genome.
#' @param win window over which to observe tandem coverage over each breakpoint
#' @param lab label for the plot
#' @param ... Other arguments passed to `get_bps`.
#'
#' @return The output is `Heatmap` object representing a coverage heatmap that
#' can be postprocessed with `smoothHeatmap(output)`, then piped into
#' `plotHeatmapList`.
#'
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps CoverageHeatmap

feature_coverage <- function(gr_ob, feat, win, lab, ...){

  # reference bps gr object
  suppressWarnings(
    ref_bps <- get_bps(gr_ob, ...) + win/2
  )

  # Remove windows that reach boundaries.
  ref_trim <- ref_bps == trim(ref_bps)
  ref_bps_trim <- ref_bps[ref_trim]

  # heatmap
  CoverageHeatmap(windows = ref_bps_trim, track = feat, coords = c(-win/2, win/2), label = lab)
}
