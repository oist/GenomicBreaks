dblExamples <- list(exampleDoubleInversion1, exampleDoubleInversion2,
                    exampleDoubleInversion3, exampleDoubleInversion4)

# Put the flanks (ranges 1 and 5) on the opposite strand, as if the query
# genome around the central triplet was read on its other strand.  The central
# triplet is unchanged, but it now requires three inversions.
flipFlanks <- function(gb) {
  strand(gb)[c(1, 5)] <- ifelse(strand(gb)[c(1, 5)] == "+", "-", "+")
  gb$query[c(1, 5)] <- gb$query[c(5, 1)]
  gb
}

test_that("each double inversion example is flagged on its left flank", {
  for (gb in dblExamples)
    expect_equal(as.vector(flagDoubleInversions(gb)$Dbl),
                 c(TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("filterDoubleInversions() returns the central triplet", {
  for (gb in dblExamples)
    expect_equal(names(filterDoubleInversions(gb)), c("2", "3", "4"))
})

test_that("removeDoubleInversions() leaves the colinear flanks", {
  for (gb in dblExamples) {
    rm <- removeDoubleInversions(gb)
    expect_length(rm, 1)
    expect_equal(as.character(granges(rm), ignore.strand = TRUE), "chrA:100-590")
    expect_equal(as.character(rm$query), "chrB:100-590")
  }
})

test_that("flagAll() sets the Dbl flag on the left flank only", {
  for (gb in dblExamples)
    expect_equal(flagAll(gb)$flag, c("Dbl", NA, NA, NA, NA))
})

test_that("flanks on the opposite strand are not flagged", {
  for (gb in dblExamples)
    expect_false(any(flagDoubleInversions(flipFlanks(gb))$Dbl))
})

test_that("flanks that are not colinear with the triplet are not flagged", {
  for (gb in dblExamples) {
    q <- gb$query
    seqlevels(q) <- c(seqlevels(q), "chrD", "chrE")
    seqnames(q)[c(1, 5)] <- c("chrD", "chrE")
    gb$query <- q
    expect_false(any(flagDoubleInversions(gb)$Dbl))
  }
})

test_that("ranges on different target sequences are not flagged", {
  for (gb in dblExamples) {
    expect_false(any(flagDoubleInversions(moveToOtherTargetSeq(gb, 1))$Dbl))
    expect_false(any(flagDoubleInversions(moveToOtherTargetSeq(gb, 5))$Dbl))
  }
})

test_that("the central triplet alone is not flagged", {
  for (gb in dblExamples)
    expect_false(any(flagDoubleInversions(gb[2:4])$Dbl))
})

test_that("short objects return a flag of the right length", {
  for (n in 0:4)
    expect_length(flagDoubleInversions(exampleDoubleInversion1[seq_len(n)])$Dbl, n)
})

test_that("other structural variant examples are not flagged", {
  for (gb in list(exampleColinear5, exampleInversion5uncollapsed,
                  exampleNestedInversions, exampleClownInversion))
    expect_false(any(flagDoubleInversions(gb)$Dbl))
})
