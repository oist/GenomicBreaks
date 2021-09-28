#' Coalesce Pairwise Alignments
#'
#' This algorithm take in a genome-to-genome alignment, represented as a
#' collection of intervals in a query genome paired with intervals in a target
#' genome.  It reduces the number of pairs by coalescing pairs that are
#' within close proximity on the same strand (user determined).
#'
#' @note Fragmented alignments arising from incorrect basecalls, misassembly or
#' misalignments can cause us to infer artificial breakpoints
#'
#' Internally, `coalesce_contigs` uses `[flagColinearAlignments()]`  See the
#' examples for details.
#'
#' @param gb [`GBreaks`] object of the pairwise alignment.
#' @param tol Unaligned region of width lesser than or equal to `tol` in both
#'        the reference and query case will be bridged in coalescing.
#' @param minwidth Remove the intervals whose width smaller than this value.
#'
#' @return Returns a new `GBreaks` object with a reduced number of alignments
#' fragments due to coalescion.  The returned object is sorted ignoring strand.
#' For convenience during analysis sessions, its _score_ is set to the width of
#' the ranges on the _target_ genome.
#'
#' @family Colinearity functions
#' @family modifier functions
#'
#' @author Charlotte West
#' @author Charles Plessy
#'
#' @examples
#' flagColinearAlignments(exampleColinear3, details = TRUE)
#' coalesce_contigs(exampleColinear3)
#'
#' # Target range [1] precedes target range [2]
#' precede(exampleColinear3)
#' # Query range [1] precedes query range [2]
#' precede(exampleColinear3$query)
#'
#' # Ranges on the minus strand
#' gb2 <- exampleColinear3 |> reverse() |> sort(ignore.strand = TRUE)
#' flagColinearAlignments(gb2, details = TRUE)
#' coalesce_contigs(gb2)
#'
#' # Target range [1] follows target range [2]
#' follow(gb2)
#' # Or, ignoring strand, target range [1] precedes target range [2]
#' precede(gb2, ignore.strand = TRUE)
#' # Query range [1] follows query range [2]
#' follow(gb2$query)
#' coalesce_contigs(gb2)
#'
#' # Ranges that should not coalesce because they are not
#' # ordered properly
#' flagColinearAlignments(exampleNotColinear, details = TRUE)
#' coalesce_contigs(exampleNotColinear)
#'
#' # Coalescing strandless objects
#' gb3 <- exampleColinear3
#' strand(gb3) <- "*"
#' flagColinearAlignments(gb3, details = TRUE)
#' coalesce_contigs(gb3)
#'
#' @export
#' @importFrom GenomicRanges GRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom stats na.omit
#' @include dist2next.R

coalesce_contigs <- function(gb, tol = Inf, minwidth = 0) {

  # Drop blocks that are narrower than `drop`

  gb <- gb[width(gb) >= minwidth]
  gb <- gb[width(gb$query) >= minwidth]

  if (length(gb) == 0)
    return(gb)

  # The rest of the algorithm assumes that the reference ranges are sorted
  gb <- sort(gb, ignore.strand = TRUE)

  # Check colinearity of query ranges
  gb <- flagColinearAlignments(gb, tol = tol, minwidth = minwidth, details = TRUE)

  # apply extension to intersected zone (applying just to end points) (ref only)
  gb$r_add <- gb$tdist
  gb[gb$colinear != TRUE]$r_add <- 0

  #######################################################################

  # define new GRanges object for output
  gr_ext <- gb
  q_ext <- gb$query
  strand(q_ext) <- strand(gb) # Guard against merging blocks on opposite strands

  end(gr_ext) <- end(gr_ext) + gb$r_add

  # reduce, concatenate, and restore original order
  reduceAndSort <- function (gr) {
    gr <- reduce(gr, min.gapwidth = 0, with.revmap = TRUE)
    gr$order <- order(unlist(lapply(gr$revmap, head, 1)))
    gr <- gr[gr$order]
    granges(gr)
  }

  gr_red <- reduceAndSort(gr_ext)

  # apply extension to intersected zone (applying just to end points) (query only)
  gb$q_add <- gb$qdist
  gb[gb$colinear != TRUE]$q_add <- 0

  end(  q_ext[strand(gb) == "+"]) <- end(  q_ext[strand(gb) == "+"]) + gb[strand(gb) == "+"]$q_add
  end(  q_ext[strand(gb) == "*"]) <- end(  q_ext[strand(gb) == "*"]) + gb[strand(gb) == "*"]$q_add
  start(q_ext[strand(gb) == "-"]) <- start(q_ext[strand(gb) == "-"]) - gb[strand(gb) == "-"]$q_add

  # reduce and concatenate
  gr_red$query <- reduceAndSort(q_ext)
  strand(gr_red$query) <- "*" # Restore un-strandedness

  score(gr_red) <- width(gr_red)

  GBreaks(sort(gr_red, ignore.strand = TRUE))
}
