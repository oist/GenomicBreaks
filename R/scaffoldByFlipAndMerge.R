#' Guided scaffolding
#'
#' @param gr A [GenomicRanges::GRanges] object.
#' @param guide A named list of data frames of contig names and orientations.
#' @param drop Drop the contigs that were not included in the guide.
#'
#' @family family scaffolding functions
#'
#' @examples
#' (gr <- exampleTranslocation |> swap())
#' (g <- list(chrBC = data.frame(contig = c("chrB", "chrC"), orientation = c(1,-1))))
#' scaffoldByFlipAndMerge(g, gr)
#'
#' @importFrom dplyr bind_rows filter pull
#'
#' @export

scaffoldByFlipAndMerge <- function(gr, guide, drop = FALSE) {
  # First, flip what should be flipped
  which.flip <- guide |>
    dplyr::bind_rows(.id='chromosome') |>
    dplyr::filter(orientation == -1) |>
    dplyr::pull(contig)
  gr <- flipContigs(gr, which.flip)

  # Then, merge what should be merged
  for (newname in names(guide)) {
    gr <- mergeSeqLevels(gr, guide[[newname]]$contig, newname)
  }
  if (isTRUE(drop)) {
    gr <- gr[seqnames(gr) %in% names(guide)]
    seqlevels(gr) <- seqlevelsInUse(gr)
  }
  sort(gr, ignore.strand = T)
}

flipContigs <- function(gr, flip_which) {
  flip_which <- intersect(flip_which, seqlevelsInUse(gr))
  if (length(flip_which) == 0) return(gr)
  grl <- split(gr, seqnames(gr))
  for (contig in flip_which)
    grl[[contig]] <- reverse(grl[[contig]])
  unlist(grl) |> unname()
}
