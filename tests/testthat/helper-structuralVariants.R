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

# Brute-force implementation of the definition of flagTranslocations(), for
# a sorted object.  Returns the expected `tra` and `trs` flags.
#  - tra: U A V consecutive on one target sequence, U and V on the same
#    strand, and V immediately after U on the query, in the direction of the
#    strand.
#  - trs: X Y consecutive on one target sequence and on the same strand, with
#    exactly one range between them on the query, in the direction of the
#    strand.
bruteTranslocations <- function(gb) {
  n  <- length(gb)
  ts <- as.character(seqnames(gb))
  qs <- as.character(seqnames(gb$query))
  st <- as.character(strand(gb))
  qrank <- ave(start(gb$query), qs, FUN = rank)
  step  <- ifelse(st == "-", -1, 1)
  tra <- trs <- logical(n)
  for (i in seq_len(n)) {
    if (i + 2 <= n) {
      U <- i; V <- i + 2
      if (ts[U] == ts[V] && st[U] == st[V] && qs[U] == qs[V] &&
          qrank[V] == qrank[U] + step[U]) tra[U] <- TRUE
    }
    if (i + 1 <= n) {
      X <- i; Y <- i + 1
      if (ts[X] == ts[Y] && st[X] == st[Y] && qs[X] == qs[Y] &&
          qrank[Y] == qrank[X] + 2 * step[X]) trs[X] <- TRUE
    }
  }
  list(tra = tra, trs = trs)
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
