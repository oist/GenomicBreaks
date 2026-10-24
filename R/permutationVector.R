#' Permutation vector
#'
#' Represent _GenomicBreaks_ objects as a _permutation vector_.
#'
#' In a `GBreaks` object properly sorted by its _target_ ranges, the _query_
#' ranges …
#'
#' @param gb A [`GBreaks`] object.
#'
#' @returns A vector of non-zero integers representing the position of the
#' genomic content in the `GBreaks` object.
#'
#' @export

permutationVector <- function(gb) {
  # Transform query ranges into simple indices
  o <- order(order(gb$query))
  o * ifelse(strand(gb) == '-', -1L, 1L)
}

#' Construct a `GBreaks` object from a permutation vector
#'
#' Inverse of [`permutationVector()`]: builds a [`GBreaks`] object of evenly
#' spaced blocks from a signed permutation.
#'
#' Block `i` of the _target_ genome is aligned to position `abs(v[i])` of the
#' _query_ genome, on the minus strand if `v[i]` is negative.  The blocks are
#' placed every `spacing` bases on each genome, starting at position `spacing`,
#' and the genomes are cut into sequences of `tsplit` and `qsplit` blocks,
#' named `chrT1`, `chrT2`, etc. and `chrQ1`, `chrQ2`, etc.  Sequence lengths are set,
#' so that the object can be passed to [`reverse()`].
#'
#' @param v A signed permutation of `1:n`, as returned by [`permutationVector()`].
#' @param tsplit,qsplit Number of blocks in each sequence of the _target_ and
#'        _query_ genomes, in order.  They must sum to `n`.
#' @param width Width of the blocks.
#' @param spacing Distance between the starts of consecutive blocks.
#'
#' @returns A `GBreaks` object of length `n`, sorted, such that
#' `permutationVector()` returns `v`.
#'
#' @seealso [`simulateRearrangements()`]
#'
#' @examples
#' # The double inversion of exampleDoubleInversion1
#' permutationVector(exampleDoubleInversion1)
#' permutationVectorToGBreaks(c(1, 4, -2, -3, 5))
#' permutationVectorToGBreaks(c(1, 4, -2, -3, 5)) |> flagDoubleInversions()
#'
#' # On two target sequences
#' permutationVectorToGBreaks(c(1, 4, -2, -3, 5), tsplit = c(2, 3))
#'
#' @export

permutationVectorToGBreaks <- function(v, tsplit = length(v), qsplit = length(v),
                                       width = 90, spacing = 100) {
  n <- length(v)
  if (!isTRUE(all.equal(sort(abs(v)), seq_len(n))))
    stop("`v` must be a signed permutation of 1:", n, ".")
  if (sum(tsplit) != n || sum(qsplit) != n || any(tsplit < 1) || any(qsplit < 1))
    stop("`tsplit` and `qsplit` must be positive and sum to ", n, ".")
  if (width >= spacing) stop("`width` must be smaller than `spacing`.")
  # Sequence name and rank of each position in a genome cut in `split` blocks
  seqOf  <- function(split, prefix) factor(rep(paste0(prefix, seq_along(split)), split),
                                           levels = paste0(prefix, seq_along(split)))
  rankOf <- function(split) sequence(split)
  seqInfo <- function(split, prefix)
    Seqinfo(paste0(prefix, seq_along(split)), (split + 1) * spacing)
  qpos <- abs(v)
  target <- GRanges(seqOf(tsplit, "chrT"),
                    IRanges(rankOf(tsplit) * spacing, width = width),
                    strand = ifelse(v < 0, "-", "+"),
                    seqinfo = seqInfo(tsplit, "chrT"))
  query  <- GRanges(seqOf(qsplit, "chrQ")[qpos],
                    IRanges(rankOf(qsplit)[qpos] * spacing, width = width),
                    seqinfo = seqInfo(qsplit, "chrQ"))
  GBreaks(target = target, query = query)
}
