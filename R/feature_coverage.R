#' Feature coverage heatmap
#'
#' Wraps the [`heatmaps::CoverageHeatmap`] function to produce a heatmap centred
#' on the boundaries of _genomic ranges_.  Assuming that these ranges are a
#' [`GBreaks`] object, then the boundaries approximate breakpoints.
#'
#' @param gr A [`GRanges`] object.
#' @param feat A `GRanges` object representing the feature of interest.
#' @param window Window over which to observe feature coverage.
#' @param label Label for the plot.
#' @param ... Other arguments passed to [`get_bps`] in order to select the
#'        boundaries, their order and their orientation.
#'
#' @return Returns a [`Heatmap`] object that can be piped into
#' [`heatmaps::smoothHeatmap`] and then [`heatmaps::plotHeatmapList`] or
#' [`heatmaps::plotHeatmapMeta`].
#'
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps CoverageHeatmap
#'
#' @family plot functions
#' @family heatmap functions
#'
#' @export


feature_coverage <- function(gr, feat, window, label, ...) {

  # reference bps gr object
  suppressWarnings(
    ref_bps <- get_bps(gr, ...) + window/2
  )

  # Remove windows that reach boundaries.
  ref_trim <- ref_bps == trim(ref_bps)
  ref_bps_trim <- ref_bps[ref_trim]

  # heatmap
  CoverageHeatmap(windows = ref_bps_trim, track = feat, coords = c(-window/2, window/2), label = label)
}
