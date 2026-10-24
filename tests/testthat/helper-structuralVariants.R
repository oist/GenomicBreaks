# Simulations stored in the package (see ?simulateRearrangements), rebuilt as
# GBreaks objects.
loadSimulations <- function(n = Inf) {
  f <- system.file("extdata/simulatedRearrangements.rds", package = "GenomicBreaks")
  # Fallback when the tests are run on sourced files, without the package
  if (f == "") f <- test_path("../../inst/extdata/simulatedRearrangements.rds")
  sims <- readRDS(f)
  sims <- head(sims, n)
  lapply(sims, \(s) permutationVectorToGBreaks(s$v, s$tsplit, s$qsplit))
}

# Brute-force implementation of the definition of flagTranslocations(), on
# one side.  The triplet X A Y is consecutive in the "main" genome, on one
# sequence, X and Y are on the same strand, and Y immediately follows X on the
# "other" genome in the direction of the strand.  Returns, for each range in
# the input order, whether it is the X of such a triplet.
bruteTranslocationSide <- function(mainSeq, mainStart, otherSeq, otherStart, strand) {
  n <- length(mainSeq)
  o <- order(mainSeq, mainStart)
  otherRank <- ave(otherStart, otherSeq, FUN = rank)
  flag <- logical(n)
  for (p in seq_len(max(n - 2, 0))) {
    X <- o[p]; A <- o[p + 1]; Y <- o[p + 2]
    if (!(mainSeq[X] == mainSeq[A] && mainSeq[A] == mainSeq[Y])) next
    if (strand[X] != strand[Y]) next
    if (otherSeq[X] != otherSeq[Y]) next
    step <- if (strand[X] == "+") 1 else -1
    if (otherRank[Y] == otherRank[X] + step) flag[X] <- TRUE
  }
  flag
}

# Expected `tra` and `query$tra` flags of flagTranslocations(gb, both = TRUE)
bruteTranslocations <- function(gb) {
  ts <- as.character(seqnames(gb));       tst <- start(gb)
  qs <- as.character(seqnames(gb$query)); qst <- start(gb$query)
  st <- as.character(strand(gb))
  list(tra  = bruteTranslocationSide(ts, tst, qs, qst, st),
       qtra = bruteTranslocationSide(qs, qst, ts, tst, st))
}

# Move the first or the last range of a GBreaks object to a new target
# sequence, ordered before or after the existing ones so that the object stays
# sorted.  The query ranges are unchanged.
moveToOtherTargetSeq <- function(gb, i) {
  stopifnot(i %in% c(1, length(gb)))
  if (i == 1) {
    seqlevels(gb) <- c("chr0", seqlevels(gb))
    seqnames(gb)[i] <- "chr0"
  } else {
    seqlevels(gb) <- c(seqlevels(gb), "chrZ")
    seqnames(gb)[i] <- "chrZ"
  }
  gb
}
