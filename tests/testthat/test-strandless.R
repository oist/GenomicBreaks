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

# Translocations: unstranded ranges are tested in both directions (experimental).

test_that("translocations in unstranded objects follow their definition", {
  for (gb in loadSimulations(60)) {
    gb <- unstrand(gb)
    f <- flagTranslocations(gb)
    b <- bruteTranslocations(gb)
    expect_equal(as.vector(f$tra), b$tra)
    expect_equal(as.vector(f$trs), b$trs)
    # Moved ranges: after each tra flag, and between the two ranges of each
    # site on the query genome, in either direction.
    qs <- as.character(seqnames(gb$query))
    qrank <- ave(start(gb$query), qs, FUN = rank)
    sitePay <- sapply(which(b$trs), \(X) {
      k <- (qrank[X + 1] - qrank[X]) / 2
      which(qs == qs[X] & qrank == qrank[X] + k)
    })
    expected <- sort(unique(c(which(b$tra) + 1, unlist(sitePay))))
    expect_equal(as.integer(names(filterTranslocations(gb))), expected)
  }
})

test_that("translocations in unstranded objects do not depend on the orientation of the query", {
  for (gb in loadSimulations(40)) {
    rev <- reverse(gb, query = TRUE) |> sort(ignore.strand = TRUE)
    f1 <- flagTranslocations(unstrand(gb))
    f2 <- flagTranslocations(unstrand(rev))
    expect_equal(as.vector(f1$tra), as.vector(f2$tra))
    expect_equal(as.vector(f1$trs), as.vector(f2$trs))
  }
})

test_that("documented example of the ambiguity of unstranded moves", {
  # Target X Y U A V, query X A Y U V
  gb <- unstrand(permutationVectorToGBreaks(c(1, 3, 4, 2, 5)))
  expect_equal(names(filterTranslocations(gb)), c("2", "3", "4"))
  # flagAll() gives colinearity the precedence, as in the oriented object.
  expect_equal(flagAll(gb)$flag, c(NA, "Col", "Tra", NA, NA))
  expect_equal(flagAll(permutationVectorToGBreaks(c(1, 3, 4, 2, 5)))$flag,
               c(NA, "Col", "Tra", NA, NA))
})

test_that("flagAll() does not fail on unstranded simulations", {
  for (gb in loadSimulations(20)) {
    expect_no_error(flagAll(unstrand(gb)))
  }
})
