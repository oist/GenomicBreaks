#' Simulate rearranged genomes
#'
#' Simulates a pair of genomes made of `n` aligned blocks: the _target_ genome
#' keeps the blocks in their ancestral order, and the _query_ genome is derived
#' from it by random inversions and transpositions.  Both genomes are then cut
#' into sequences at random positions.
#'
#' Each event picks its blocks at random positions of the _query_ genome as it
#' stands after the previous events, before it is cut into sequences:
#'
#'  * `inversion`: reverses a segment of one block or more and flips its
#'    strand.
#'  * `transposition`: moves a segment of `1` to `maxMove` blocks to another
#'    position.
#'  * `inverted transposition`: same, and inverts the segment.
#'
#' Events are not checked for parsimony: a later event can undo, or extend, the
#' effect of an earlier one.  The random cuts of each genome into sequences are
#' not recorded as events.
#'
#' @param n Number of blocks.
#' @param nEvents Number of events.
#' @param types Types of events to sample from (see Details).
#' @param maxMove Maximum number of blocks moved by a transposition.
#' @param nTargetSeqs,nQuerySeqs Number of sequences in each genome.
#' @param ... Passed to [`permutationVectorToGBreaks()`] (`width`, `spacing`).
#'
#' @returns A sorted [`GBreaks`] object of length `n`.  Its `metadata()` list
#' contains an `events` data frame with one row per event, in order: `type`;
#' `segment`, the blocks affected by the event, in _query_ order before the
#' event, with a minus sign for the blocks on the minus strand; and `after`,
#' for transpositions, the block preceding the moved segment in the _query_
#' genome after the event (`0` when the segment was moved to the start).
#' Blocks are identified by their position in the _target_ genome, which is
#' the ancestral order.
#'
#' @seealso [`permutationVectorToGBreaks()`] and [`permutationVector()`].  A set
#' of 200 simulations is provided in the `simulatedRearrangements.rds` file of
#' the package (see Examples).
#'
#' @examples
#' set.seed(1)
#' gb <- simulateRearrangements(n = 10, nEvents = 2)
#' gb
#' metadata(gb)$events
#' permutationVector(gb)
#' flagAll(gb)
#'
#' # Pre-computed simulations, stored as permutation vectors
#' sims <- readRDS(system.file("extdata/simulatedRearrangements.rds",
#'                             package = "GenomicBreaks"))
#' sims[[1]]
#' with(sims[[1]], permutationVectorToGBreaks(v, tsplit, qsplit))
#'
#' @importFrom S4Vectors metadata metadata<-
#' @export

simulateRearrangements <- function(n = 14, nEvents = 3,
                                   types = c("inversion", "transposition", "inverted transposition"),
                                   maxMove = 2, nTargetSeqs = 1, nQuerySeqs = 1, ...) {
  types <- match.arg(types, several.ok = TRUE)
  if (n < maxMove + 2) stop("`n` must be at least `maxMove + 2`.")
  if (nTargetSeqs > n || nQuerySeqs > n) stop("There can not be more sequences than blocks.")
  # Signed ancestral block identifiers, in query order
  q <- seq_len(n)
  fmt <- function(x) paste(x, collapse = " ")
  events <- data.frame(type = character(nEvents), segment = character(nEvents),
                       after = rep(NA_integer_, nEvents))
  for (k in seq_len(nEvents)) {
    type <- if (length(types) == 1) types else sample(types, 1)
    if (type == "inversion") {
      ij  <- sort(sample(n, 2, replace = TRUE))
      pos <- ij[1]:ij[2]
      events$segment[k] <- fmt(q[pos])
      q[pos] <- -rev(q[pos])
    } else {
      len  <- if (maxMove == 1) 1 else sample(maxMove, 1)
      i    <- sample(n - len + 1, 1)
      pos  <- i:(i + len - 1)
      seg  <- q[pos]
      rest <- q[-pos]
      events$segment[k] <- fmt(seg)
      if (type == "inverted transposition") seg <- -rev(seg)
      # Exclude the original position, where the event would do nothing (or
      # would be an inversion).
      dest <- setdiff(0:length(rest), i - 1)
      j    <- if (length(dest) == 1) dest else sample(dest, 1)
      q    <- append(rest, seg, after = j)
      events$after[k] <- if (j == 0) 0L else as.integer(rest[j])
    }
    events$type[k] <- type
  }
  # Cut a genome of n blocks into s sequences at random positions
  randomSplit <- function(s) diff(c(0, sort(sample(n - 1, s - 1)), n))
  # Signed query position of each target block
  v <- integer(n)
  v[abs(q)] <- as.integer(sign(q)) * seq_len(n)
  gb <- permutationVectorToGBreaks(v, tsplit = randomSplit(nTargetSeqs),
                                   qsplit = randomSplit(nQuerySeqs), ...)
  metadata(gb)$events <- events
  gb
}
