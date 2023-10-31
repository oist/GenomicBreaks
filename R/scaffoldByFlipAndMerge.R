#' Guided scaffolding
#'
#' @param gr A [GenomicRanges::GRanges] object.
#' @param guide A named list of data frames of contig names and orientations.
#' @param drop Drop the contigs that were not included in the guide.  The new
#'        sequence levels will be in the same order as in the `guide`.
#'
#' @family scaffolding functions
#'
#' @author Michael Mansfield
#' @author Charles Plessy
#'
#' @examples
#' (gr <- GRanges(c("chrB:100-200:+", "chrC:201-300:+",
#'                  "chrB:301-400:+", "chrD:100-200:+")) |> forceSeqLengths())
#' (g <- list(
#'   chrD = data.frame(contig = "chrD", orientation = 1),
#'   chrBC = data.frame(contig = c("chrB", "chrC"), orientation = c(1,-1))))
#'
#' scaffoldByFlipAndMerge(gr, g, drop = TRUE)
#'
#' @export

scaffoldByFlipAndMerge <- function(gr, guide, drop = FALSE) {
  # First, flip what should be flipped
  which.flip <- do.call(rbind, guide)
  which.flip <- which.flip[which.flip$orientation == -1, "contig", drop = TRUE]
  gr <- flipContigs(gr, which.flip)

  # Then, merge what should be merged
  for (newname in names(guide)) {
    gr <- mergeSeqLevels(gr, guide[[newname]]$contig, newname)
  }
  if (isTRUE(drop)) {
    gr <- gr[seqnames(gr) %in% names(guide)]
    seqlevels(gr) <- names(guide)
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
