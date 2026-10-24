test_that("the gap of the bundled example", {
  lg <- leftInversionGaps(exampleInversion)
  expect_equal(as.character(lg), "chrA:191-209")
  expect_equal(lg$flagIndex, 1)
})

test_that("one gap per inversion, including zero-width gaps", {
  # The left flank and the central range of the second inversion are adjacent
  # on the target genome.
  gb <- GBreaks(target = GRanges(c("chrA:100-190:+", "chrA:210-291:-", "chrA:301-400:+",
                                   "chrA:500-590:+", "chrA:591-690:-", "chrA:701-800:+")),
                query  = GRanges(c("chrB:100-190",   "chrB:210-291",   "chrB:301-400",
                                   "chrB:500-590",   "chrB:600-690",   "chrB:701-800")))
  lg <- leftInversionGaps(gb)
  expect_equal(lg$flagIndex, c(1, 4))
  expect_equal(start(lg), c(191, 591))
  expect_equal(width(lg), c(19, 0))
})

test_that("detect = FALSE reuses the existing flags", {
  expect_length(leftInversionGaps(flagInversions(exampleInversion, tol = 18), detect = FALSE), 0)
  expect_length(leftInversionGaps(flagInversions(exampleInversion, tol = 19), detect = FALSE), 1)
  expect_length(leftInversionGaps(exampleInversion, detect = FALSE), 0)
})

test_that("objects without inversions return no gap", {
  expect_length(leftInversionGaps(exampleColinear3), 0)
  expect_length(leftInversionGaps(exampleColinear3[0]), 0)
})
