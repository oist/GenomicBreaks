#' Flag colinear alignments
#'
#' Flags alignments that are colinear with the next one in sequence order.
#' The flag is added to the first alignment.  Colinearity is defined by the fact
#' that the next alignment on the _target_ genome is the same as the next
#' alignment on the _query_ genome, with strand information being taken into
#' account.
#'
#' Internally, `flagColinearAlignments()` uses the [`GenomicRanges::precede`]
#' and [`GenomicRanges::follow`] functions functions to determine what is the
#' _next_.  For a given range, these functions return the index position of the
#' range it precedes or follows, or `NA` as the first range follows nothing and
#' the last range precedes nothing.  See the examples for details.
#'
#' @note The flags are only valid as long as the Genomic Breaks object is not
#' modified by removing alignments or sorting them in a different order.
#'
#' @param gb [`GBreaks`] object of the pairwise alignment.
#' @param tol Unaligned regions larger than this _tolerance_ threshold will
#'        interrupt colinearity.
#' @param minwidth Remove the intervals whose width smaller than this value.
#' @param details Returns more metadata columns if `TRUE`.
#'
#' @note Pay attention that if the `mindwidth` option is passed, some intervals
#' are discarded from the returned object.  This parameter might be removed in
#' the future if too confusing or useless.
#'
#' @return Returns a modified `GBreaks` object with a new `colinear` metadata
#' column indicating if an alignment is colinear with the next one.  If the
#' `details` option is set to `TRUE`, it will also output the relative position
#' of the following and previous alignments on the _target_ and _query_ genomes
#' (`tfoll`, `tprev`, `qfoll`, `qprev`), a partial colinearity flag for each
#' genome (`t_col` and `q_col`), and the distance to the next alignment on
#' each genome (`tdist` and `qdist`).
#'
#' @examples
#' # Ranges on the plus strand that should coalesce
#' gb1       <- GenomicRanges::GRanges(c(A="Ref:100-200:+", B="Ref:400-500:+"))
#' gb1$query <- GenomicRanges::GRanges(c(A="Que:100-200",   B="Que:400-500"))
#' gb1
#' flagColinearAlignments(gb1)
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
#' flagColinearAlignments(gb2)
#'
#' # Ranges on the minus strand that should not coalesce because they are not
#' # ordered properly
#' gb3       <- GenomicRanges::GRanges(c("Ref:100-200:-", "Ref:400-500:-"))
#' gb3$query <- GenomicRanges::GRanges(c("Que:100-200",   "Que:400-500"))
#' # Reference range [1] follows reference range [2]
#' GenomicRanges::follow(gb3)
#' # Query range [1] follows query range [2]
#' GenomicRanges::precede(gb3$query)
#' flagColinearAlignments(gb3)
#'
#' # Ranges on the plus strand that should not coalesce because they are not
#' # ordered properly
#' gb4       <- GenomicRanges::GRanges(c("Ref:100-200:+", "Ref:400-500:+", "Ref:600-700:+"))
#' gb4$query <- GenomicRanges::GRanges(c("Que:1100-1200:+", "Que:1700-1800:+", "Que:1500-1600:+"))
#' flagColinearAlignments(gb4)
#'
#' @family Flagging functions
#' @family Colinearity functions
#'
#' @author Charlotte West
#' @author Charles Plessy
#'
#' @export
#' @importFrom GenomicRanges GRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @importFrom stats na.omit
#' @include dist2next.R

flagColinearAlignments <- function(gb, tol = Inf, minwidth = 0, details = FALSE) {
  # Drop blocks that are narrower than `minwidth`
  gb <- gb[width(gb) >= minwidth]
  gb <- gb[width(gb$query) >= minwidth]

  # Handle empty objects
  if (length(gb) == 0) return(gb)

  # Sort the ranges following the target genome.  The argument
  # ignore.strand = TRUE is set so that alignments on the opposite strand can
  # interrupt colinearity.
  gb <- sort(gb, ignore.strand = TRUE)

  # Check colinearity of query ranges

  # The names of the precede() and follow() functions are a bit confusing;
  # check the help page if needed.

  # Position of the next block minus position of the current block equals
  # to 1 when they are colinear.
  gb$tfoll <- precede(gb)       - seq_along(gb)
  gb$qfoll <- precede(gb$query) - seq_along(gb)

  # Position of the previous block minus position of the current block equals
  # to 1 when they are anti-colinear.
  gb$tprev <- follow( gb)       - seq_along(gb)
  gb$qprev <- follow( gb$query) - seq_along(gb)

  # When the reference strand is "+", we want the query blocks to be colinear
  # and when the reference is "-" we want them to be anti-colinear.
  gb$t_col <- ( strand(gb) == "+" & gb$tfoll == 1 ) |
              ( strand(gb) == "-" & gb$tprev == 1 )

  gb$q_col <- ( strand(gb) == "+" & gb$qfoll == 1 ) |
              ( strand(gb) == "-" & gb$qprev == 1 )

  # Calculate distance and check if it is within tolerance threshold
  # See 18939aa98b2db83a2360d5a0e92966581de4799e about tol + 1
  gb <- dist2next(gb)
  gb$colinear <- gb$tdist < tol + 1 &
                 gb$qdist < tol + 1 &
                 gb$q_col
  gb$colinear[is.na(gb$colinear)] <- FALSE

  if (! isTRUE(details)) {
    gb$tfoll <- gb$qfoll <- gb$tprev <- gb$qprev <- NULL
    gb$t_col <- gb$q_col <- gb$tdist <- gb$qdist <- NULL
  }

  gb
}
