test_that("permutationVectorToGBreaks() is the inverse of permutationVector()", {
  set.seed(1)
  for (k in 1:20) {
    n <- sample(2:12, 1)
    v <- sample(n) * sample(c(-1L, 1L), n, replace = TRUE)
    tsplit <- diff(c(0, sort(sample(n - 1, sample(0:min(2, n - 1), 1))), n))
    qsplit <- diff(c(0, sort(sample(n - 1, sample(0:min(2, n - 1), 1))), n))
    gb <- permutationVectorToGBreaks(v, tsplit, qsplit)
    expect_true(isSorted(gb))
    expect_equal(permutationVector(gb), v)
    expect_equal(as.vector(table(seqnames(gb))), tsplit)
    expect_equal(as.vector(table(seqnames(gb$query))), qsplit)
  }
})

test_that("permutationVectorToGBreaks() builds the bundled double inversion", {
  gb <- permutationVectorToGBreaks(c(1, 4, -2, -3, 5))
  expect_equal(permutationVector(gb), permutationVector(exampleDoubleInversion1))
  expect_equal(as.vector(flagDoubleInversions(gb)$Dbl), as.vector(flagDoubleInversions(exampleDoubleInversion1)$Dbl))
  expect_false(anyNA(seqlengths(gb)))
})

test_that("permutationVectorToGBreaks() rejects invalid input", {
  expect_error(permutationVectorToGBreaks(c(1, 1, 3)))
  expect_error(permutationVectorToGBreaks(1:3, tsplit = c(1, 1)))
  expect_error(permutationVectorToGBreaks(1:3, width = 100))
})

test_that("simulateRearrangements() returns a sorted object and its events", {
  set.seed(2)
  gb <- simulateRearrangements(n = 12, nEvents = 4, nTargetSeqs = 2, nQuerySeqs = 3)
  expect_s4_class(gb, "GBreaks")
  expect_length(gb, 12)
  expect_true(isSorted(gb))
  expect_equal(nlevels(seqnames(gb)), 2)
  expect_equal(nlevels(seqnames(gb$query)), 3)
  ev <- metadata(gb)$events
  expect_equal(nrow(ev), 4)
  expect_true(all(ev$type %in% c("inversion", "transposition", "inverted transposition")))
  expect_true(all(is.na(ev$after) == (ev$type == "inversion")))
})

test_that("simulateRearrangements() without events is colinear", {
  gb <- simulateRearrangements(n = 8, nEvents = 0)
  expect_equal(permutationVector(gb), 1:8)
  expect_length(coalesce_contigs(gb), 1)
  expect_equal(nrow(metadata(gb)$events), 0)
})

test_that("one inversion inside the genome is flagged as an inversion", {
  set.seed(3)
  for (k in 1:10) {
    gb <- simulateRearrangements(n = 10, nEvents = 1, types = "inversion")
    seg <- as.integer(strsplit(metadata(gb)$events$segment, " ")[[1]])
    if (1 %in% seg || 10 %in% seg) next  # An inversion at an end has no left or right flank
    expect_equal(sum(flagInversions(coalesce_contigs(gb))$inv), 1)
  }
})
