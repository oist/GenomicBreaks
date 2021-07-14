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
#' Internally, the `coalesce_contigs` function uses the `precede` and `follow`
#' functions of the `GenomicRanges` package.  For a given range, these functions
#' return the index position of the range it precedes or follows, or `NA` as the
#' first range follows nothing and the last range precedes nothing.  See the
#' examples for details.
#'
#' @param gb [`GBreaks`] object of the pairwise alignment.
#' @param tol Unaligned region of width lesser than or equal to `tol` in both
#'        the reference and query case will be bridged in coalescing.
#' @param minwidth Remove the intervals whose width smaller than this value.
#'
#' @return Returns a new `GBreaks` object with a reduced number of alignments
#' fragments due to coalescion.  The returned object is sorted ignoring strand.
#'
#' @family Colinearity functions
#'
#' @author Charlotte West
#' @author Charles Plessy
#'
#' @examples
#' # Ranges on the plus strand that should coalesce
#' gb1       <- GenomicRanges::GRanges(c(A="Ref:100-200:+", B="Ref:400-500:+"))
#' gb1$query <- GenomicRanges::GRanges(c(A="Que:100-200",   B="Que:400-500"))
#' gb1
#' coalesce_contigs(gb1)
#'
#' # Reference range [1] precedes reference range [2]
#' GenomicRanges::precede(gb1)
#' # Query range [1] precedes query range [2]
#' GenomicRanges::precede(gb1$query)
#'
#' # Ranges on the minus strand that should coalesce
#' gb2       <- GenomicRanges::GRanges(c(B="Ref:100-200:-", A="Ref:400-500:-"))
#' gb2$query <- GenomicRanges::GRanges(c(B="Que:400-500",   A="Que:100-200"))
#' gb2$qname <- names(gb2$query)
#' gb2
#' # Reference range [1] follows reference range [2]
#' GenomicRanges::follow(gb2)
#' # Or, ignoring strand, reference range [1] precedes reference range [2]
#' GenomicRanges::precede(gb2, ignore.strand = TRUE)
#' # Query range [1] follows query range [2]
#' GenomicRanges::follow(gb2$query)
#' coalesce_contigs(gb2)
#'
#' # Ranges on the minus strand that should not coalesce because they are not
#' # ordered properly
#' gb3       <- GenomicRanges::GRanges(c("Ref:100-200:-", "Ref:400-500:-"))
#' gb3$query <- GenomicRanges::GRanges(c("Que:100-200",   "Que:400-500"))
#' # Reference range [1] follows reference range [2]
#' GenomicRanges::follow(gb3)
#' # Query range [1] follows query range [2]
#' GenomicRanges::precede(gb3$query)
#' coalesce_contigs(gb3)
#'
#' Ranges on the plus strand that should not coalesce because they are not
#' ordered properly
#' gb4       <- GenomicRanges::GRanges(c("Ref:100-200:+", "Ref:400-500:+", "Ref:600-700:+"))
#' gb4$query <- GenomicRanges::GRanges(c("Que:1100-1200:+", "Que:1700-1800:+", "Que:1500-1600:+"))
#' coalesce_contigs(gb4)
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

  # The names of the precede() and follow() functions are a bit confusing;
  # check the help page if needed.

  # Position of the next block minus position of the current block equals
  # to 1 when they are colinear.
  gb$qnext <- precede(gb$query) - seq_along(gb$query)

  # Position of the previous block minus position of the current block equals
  # to 1 when they are anti-colinear.
  gb$qprev <- follow( gb$query) - seq_along(gb$query)

  # When the reference strand is "+", we want the query blocks to be colinear
  # and when the reference is "-" we want them to be anti-colinear.
  gb$q_col_with_next <- ( strand(gb) == "+" & gb$qnext == 1 ) |
                           ( strand(gb) == "-" & gb$qprev == 1 )

  gb <- dist2next(gb)

  #######################################################################

  # find intersection
  gb$con_met_total <- gb$tdist < tol + 1 &
                         gb$qdist < tol + 1 &
                         gb$q_col_with_next
  gb$con_met_total[is.na(gb$con_met_total)] <- FALSE

  # apply extension to intersected zone (applying just to end points) (ref only)
  gb$r_add <- gb$tdist
  gb[gb$con_met_total != TRUE]$r_add <- 0

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
  gb[gb$con_met_total != TRUE]$q_add <- 0

  end(  q_ext[strand(gb) == "+"]) <- end(  q_ext[strand(gb) == "+"]) + gb[strand(gb) == "+"]$q_add
  start(q_ext[strand(gb) == "-"]) <- start(q_ext[strand(gb) == "-"]) - gb[strand(gb) == "-"]$q_add

  # reduce and concatenate
  gr_red$query <- reduceAndSort(q_ext)
  strand(gr_red$query) <- "*" # Restore un-strandedness

  sort(gr_red, ignore.strand = TRUE)
}
