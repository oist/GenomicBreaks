#' Sliding window
#'
#' This function is used to subset a `GBreaks` object into a list of `GBreaks`
#' using a sliding window.
#'
#' @param gb A [`GBreaks`] object.
#' @param windowSize The size of the window.
#' @param stepSize The size of the step.
#' @param type slide the window on the reference genome or query genome (reference/query)
#' @param merged set TRUE to output a ['GBreaks'] object
#' @param cut whether cut the alignments out of the window
#'
#' @returns Returns a list of `GBreaks` objects, one for each window.
#'
#' @family Similarity indexes
#'
#' @author Zhang Kun
#'
#' @examples
#' exdata_Sac <- system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks")
#' gb <- load_genomic_breaks(exdata_Sac, BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae)
#' slided_gb <- slidingWindow(gb, windowSize = 1e5, stepSize = 5e4)
#' slided_gb[[1]]
#'
#' @export

slidingWindow <- function(gb,
                          factor = 100,
                          windowSize = NULL,
                          stepSize = NULL,
                          merged = FALSE,
                          cut = TRUE
) {

  gb <- forceSeqLengths(gb)

  if (is.null(windowSize) && is.null(stepSize)) {
    total_len <- sum(seqlengths(gb), na.rm = TRUE)
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

  tiles <- seqinfo(gb) |> as("GRanges") |> slidingWindows(windowSize, stepSize) |> unlist(use.names = FALSE)
  tiles <- tiles[width(tiles) == windowSize] # Discard the last tile if too short.

  if (is.null(tiles)) {
    return(GRangesList())
  }

  hits <- GenomicRanges::findOverlaps(gb, tiles)
  gb_hits <- gb[queryHits(hits)]
  gb_hits$overlapping_window <- tiles[subjectHits(hits)]
  gb_hits$window_id          <-       subjectHits(hits)

  if (merged) return(gb_hits)

  gb_list <- split(gb_hits, gb_hits$window_id)
  if (!cut) return(gb_list)

  trim_func <- function(gr) {
    window <- unique(gr$overlapping_window)

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
    q_strand <- as.character(strand(query))
    q_width <- width(query)
    ref_width <- width(gr)

    new_q_start <- q_start
    new_q_end   <- q_end

    is_fwd <- q_strand == "+"
    is_rev <- q_strand == "-"

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
    score(gr_ref_cut) <- round(score(gr) * ( width(gr_ref_cut) / ref_width ))
    gr_ref_cut$query <- new_query
    gr_ref_cut
  }

  lapply(gb_list, trim_func) |> GRangesList()
}
