test_that("flagTranslocations() follows its definition", {
  ntra <- ntrs <- 0
  for (gb in loadSimulations()) {
    f <- flagTranslocations(gb)
    b <- bruteTranslocations(gb)
    expect_equal(as.vector(f$tra), b$tra)
    expect_equal(as.vector(f$trs), b$trs)
    ntra <- ntra + sum(b$tra); ntrs <- ntrs + sum(b$trs)
  }
  # The test is not void
  expect_gt(ntra, 0)
  expect_gt(ntrs, 0)
})

test_that("a move along one sequence is flagged at its payload region and at its site", {
  # Target X Y U A V, query X A Y U V
  gb <- permutationVectorToGBreaks(c(1, 3, 4, 2, 5))
  f <- flagTranslocations(gb)
  expect_equal(as.vector(f$tra), c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(as.vector(f$trs), c(TRUE,  FALSE, FALSE, FALSE, FALSE))
  # Same on the minus strand: the flags stay on the first range in target order
  f <- flagTranslocations(reverse(gb, query = TRUE) |> sort(ignore.strand = TRUE))
  expect_equal(as.vector(f$tra), c(FALSE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(as.vector(f$trs), c(TRUE,  FALSE, FALSE, FALSE, FALSE))
})

test_that("a move seen in one genome only is flagged once", {
  # The moved range is alone on its query sequence, so that only the target
  # genome shows a context for it; in the swapped object, only the query
  # genome does.
  f <- flagTranslocations(exampleTranslocation)
  expect_equal(as.vector(f$tra), c(TRUE,  FALSE, FALSE))
  expect_equal(as.vector(f$trs), c(FALSE, FALSE, FALSE))
  f <- flagTranslocations(swap(exampleTranslocation, sort = TRUE))
  expect_equal(as.vector(f$tra), c(FALSE, FALSE, FALSE))
  expect_equal(as.vector(f$trs), c(TRUE,  FALSE, FALSE))
})

test_that("inversions and colinear ranges are not flagged", {
  for (gb in list(exampleColinear3, exampleInversion, exampleDoubleInversion1,
                  exampleTwinInversions)) {
    f <- flagTranslocations(gb)
    expect_false(any(f$tra))
    expect_false(any(f$trs))
  }
})

test_that("short objects return flags of the right length", {
  for (n in 0:2) {
    f <- flagTranslocations(exampleTranslocation[seq_len(n)])
    expect_length(f$tra, n)
    expect_length(f$trs, n)
  }
})

test_that("filterTranslocations() returns or removes the central range", {
  expect_equal(names(filterTranslocations(exampleTranslocation)), "2")
  expect_equal(names(filterTranslocations(exampleTranslocation, remove = TRUE)), c("1", "3"))
})

test_that("filterTranslocations(remove = TRUE) keeps everything when nothing is flagged", {
  expect_length(filterTranslocations(exampleColinear3), 0)
  kept <- filterTranslocations(exampleColinear3, remove = TRUE)
  expect_equal(names(kept), c("1", "2", "3"))
  expect_equal(granges(unname(kept)), granges(exampleColinear3))
})
