#' Breakpoint Associated Heatmaps
#'
#' This function generates a heatmap of the specified pattern, over breakpoints
#' aligned at the centre of the plot. The result can be piped into
#' [`heatmaps::smoothHeatmap`] and then [`heatmaps::plotHeatmapList`] or
#' [`heatmaps::plotHeatmapMeta`]
#'
#' @param gr `GRanges` object containing pairwise alignment
#' @param window Range over which heatmap characteristic is plotted.
#'        Breakpoints will be aligned at the center of this.
#' @param pattern Character string of desired pattern/characteristic to be
#'        plotted on heatmap.
#' @param ... Pass other arguments to [`get_bps`].
#'
#' @return Returns a [`heatmaps::Heatmap`] object of `pattern` around centred breakpoints.
#'
#' @note The `GRanges` object is expected to have a _sequence information_
#' (see [`Seqinfo::seqinfo`]) that allows the retrieval its corresponding `BSgenome`
#' object via the [`BSgenome::getBSgenome`] function.
#'
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom heatmaps PatternHeatmap
#' @importFrom Biostrings getSeq
#'
#' @author Charlotte West
#' @author Charles Plessy
#'
#' @family plot functions
#' @family heatmap functions
#'
#' @examples
#' # The plot makes no sense, but that is the best example I have at the moment.
#' exdata_Sac <- system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks")
#' gb <- load_genomic_breaks(exdata_Sac, BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae)
#' bp_heatmap(gb, 200, 'GC', dir = "left") |>
#'   heatmaps::smoothHeatmap() |> heatmaps::plotHeatmapList()
#'
#' @export

bp_heatmap <- function (gr, window, pattern, ...) {

  # Suppress warnings about overflow
  suppressWarnings(
  gr_bps <- get_bps(gr, ...) + window / 2
  )
  # Remove the out-of-bound ranges.
  gr_bps <- gr_bps[gr_bps == trim(gr_bps)]

  # Detect genome name and get BSgenome object
  genomeName <- unique(genome(gr))
  if(is.na(genomeName)) stop("No BSgenome information in the GRanges/GBreaks object.")
  genome <- BSgenome::getBSgenome(genomeName)

  # Heatmap
  hm_prep <- getSeq(genome, gr_bps)
  PatternHeatmap(hm_prep, pattern = pattern, c(-window / 2, window / 2))
}
