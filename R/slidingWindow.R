#' Sliding window
#'
#' This function is used to subset a `GBreaks` object into a list of `GBreaks`
#' using a sliding window.
#'
#' @param gb A [`GBreaks`] object.
#' @param factor Number of windows.
#' @param windowSize The size of the window.
#' @param stepSize The size of the step.
#' @param merged set TRUE to output a `GBreaks` object
#' @param cut whether cut the alignments out of the window
#'
#' @returns Returns a [`GenomicRanges::GRangesList`] containing one `GBreaks`
#' object for each window, or a single `GBreaks` object if `merged` is `TRUE`.
#'
#' @details When `cut` is `TRUE`, an alignment that crosses a window boundary
#' is cut at that boundary on the _target_ genome.  Its _query_ range is cut by
#' the same proportion, on the side given by the strand of the _target_ range,
#' and its score, if any, is scaled by that proportion.
#'
#' @section Unstranded ranges (experimental):
#' The _query_ ranges of alignments on the `*` strand are not cut.
#'
#' @family Similarity indexes
#'
#' @author Zhang Kun
#'
#' @examples
#' exampleInversion
#' slided_gb <- slidingWindow(exampleInversion, windowSize = 250, stepSize = 125)
#' slided_gb
#' sapply(slided_gb, strand_randomisation_index)
#'
#' @export

slidingWindow <- function(gb,
                          factor = 100,
                          windowSize = NULL,
                          stepSize = NULL,
                          merged = FALSE,
                          cut = TRUE
) {

  has_seqlengths <- function(gr) {
    any(!is.na(seqlengths(gr)))
  }
  assign_seqlengths_from_max <- function(gr) {
    # Ensure input is a GRanges or GBreaks object
    if (!inherits(gr, "GRanges") && !inherits(gr, "GBreaks")) {
      stop("Input must be a GRanges or GBreaks object")
    }
    # Compute maximum end per chromosome
    max_ends <- tapply(end(gr), as.character(seqnames(gr)), max, na.rm = TRUE)
    # Assign as seqlengths
    seqlengths(gr)[names(max_ends)] <- max_ends
    return(gr)
  }
  if (!has_seqlengths(gb)) {
    gb <- assign_seqlengths_from_max(gb)
  }

  if (is.null(windowSize) && is.null(stepSize)) {
    seqlens   <- seqlengths(gb)
    total_len <- sum(seqlens, na.rm = TRUE)
    # window size ≈ genome length / # windows
    windowSize <- floor(
      total_len * 2 / (factor + 1)
    )
    stepSize <- floor(windowSize / 2)
    message("Computed windowSize = ", windowSize,
            "  stepSize = ", stepSize,
            " to get ~", factor, " windows.")
  } else {
    if (!is.null(windowSize) && windowSize %% 1 != 0 || windowSize <= 0) stop("windowSize must be a positive integer")
    if (!is.null(stepSize) && stepSize %% 1 != 0 || stepSize <= 0) stop("stepSize must be a positive integer")
  }

  makeSlidingWindows <- function(seqlengths_vec, windowSize, stepSize) {
    gr_list <- lapply(names(seqlengths_vec), function(chr) {
      chr_len <- as.numeric(seqlengths_vec[[chr]])
      if (chr_len < windowSize) {
        return(GRanges(seqnames = chr, ranges = IRanges(start = 1, end = chr_len)))
      }
      starts <- seq(1, chr_len - windowSize + 1, by = stepSize)
      ends <- pmin(starts + windowSize - 1, chr_len)
      GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
    })
    do.call(c, gr_list)
  }

  seqlengths_vec <- seqlengths(gb)
  tiles <- makeSlidingWindows(seqlengths_vec, windowSize, stepSize)

  if (is.null(tiles)) {
    return(NULL)
  }

  hits <- GenomicRanges::findOverlaps(gb, tiles)
  gb_hits <- gb[queryHits(hits)]
  overlapping_windows <- tiles[subjectHits(hits)]

  mcols(gb_hits)$window_id     <- subjectHits(hits)
  mcols(gb_hits)$window_seqname <- as.character(seqnames(overlapping_windows))
  mcols(gb_hits)$window_start   <- start(overlapping_windows)
  mcols(gb_hits)$window_end     <- end(overlapping_windows)

  if (merged) return(gb_hits)

  # Cuts each range of `gr` to the parallel range of `window`.  All operations
  # are element-wise, so all the alignment-window pairs are cut in one call.
  trim_func <- function(gr, window) {
    gr_ref_cut <- GenomicRanges::pintersect(gr, window)

    ref_start_orig <- start(gr)
    ref_end_orig   <- end(gr)
    ref_start_cut  <- start(gr_ref_cut)
    ref_end_cut    <- end(gr_ref_cut)

    cut_5p <- pmax(0, ref_start_cut - ref_start_orig)
    cut_3p <- pmax(0, ref_end_orig - ref_end_cut)

    query <- gr$query
    q_start <- start(query)
    q_end   <- end(query)
    q_width <- width(query)
    ref_width <- width(gr)

    new_q_start <- q_start
    new_q_end   <- q_end

    # The orientation of an alignment is the strand of its target range.
    is_fwd <- as.character(strand(gr)) == "+"
    is_rev <- as.character(strand(gr)) == "-"

    new_q_start[is_fwd] <- q_start[is_fwd] + round(q_width[is_fwd] * cut_5p[is_fwd] / ref_width[is_fwd])
    new_q_end[is_fwd]   <- q_end[is_fwd]   - round(q_width[is_fwd] * cut_3p[is_fwd] / ref_width[is_fwd])
    new_q_start[is_rev] <- q_start[is_rev] + round(q_width[is_rev] * cut_3p[is_rev] / ref_width[is_rev])
    new_q_end[is_rev]   <- q_end[is_rev]   - round(q_width[is_rev] * cut_5p[is_rev] / ref_width[is_rev])

    new_q_start <- pmin(new_q_start, new_q_end)
    new_q_end   <- pmax(new_q_start, new_q_end)

    new_query <- GRanges(
      seqnames = seqnames(query),
      ranges = IRanges(start = new_q_start, end = new_q_end),
      strand = strand(query),
      seqinfo = seqinfo(query)
    )

    mcols(gr_ref_cut) <- mcols(gr)
    gr_ref_cut$query <- new_query
    # Scale the alignment scores, if any, to the retained fraction of the target
    if (!is.null(gr$score))
      gr_ref_cut$score <- round(gr$score * width(gr_ref_cut) / ref_width)
    gr_ref_cut
  }

  if (cut) gb_hits <- trim_func(gb_hits, overlapping_windows)

  split(gb_hits, gb_hits$window_id)
}
