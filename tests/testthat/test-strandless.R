# Inversions can not be determined without orientation.

test_that("unstranded objects have no inversions", {
  for (gb in list(exampleInversion, exampleNestedInversions, exampleDoubleInversion1,
                  exampleDoubleInversion2, exampleDoubleInversion3, exampleDoubleInversion4,
                  exampleTwinInversions)) {
    gb <- unstrand(gb)
    expect_false(any(flagInversions(gb)$inv))
    expect_false(any(flagDoubleInversions(gb)$Dbl))
    expect_false(any(flagTwinInversions(gb)$twi))
  }
})

test_that("a triplet with an unstranded range is not an inversion", {
  # Colinear, with an unstranded central range
  gb <- GBreaks(target = GRanges(c("chrA:100-190:+", "chrA:200-290:*", "chrA:300-390:+")),
                query  = GRanges(c("chrB:100-190",   "chrB:200-290",   "chrB:300-390")))
  expect_false(any(flagInversions(gb)$inv))
  # Inverted central range, with unstranded flanks
  strand(gb) <- c("*", "-", "*")
  expect_false(any(flagInversions(gb)$inv))
})

test_that("bridge regions do not create inversions and do not hide them", {
  gb <- withBridges(exampleColinear5)
  expect_true(any(strand(gb) == "*"))
  expect_false(any(flagInversions(gb)$inv))
  expect_no_error(flagAll(gb))
  # Inversion of range 3, with bridges between ranges 1 and 2, and 4 and 5.
  gb <- withBridges(permutationVectorToGBreaks(c(1, 2, -3, 4, 5)))
  expect_equal(as.character(strand(gb)), c("+", "*", "+", "-", "+", "*", "+"))
  expect_equal(which(as.vector(flagInversions(gb)$inv)), 3)
})
