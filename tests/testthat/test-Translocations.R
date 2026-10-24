test_that("flagTranslocations() follows its definition on both genomes", {
  nflags <- 0
  for (gb in loadSimulations(25)) {
    f <- flagTranslocations(gb)
    b <- bruteTranslocations(gb)
    expect_equal(as.vector(f$tra),       b$tra)
    expect_equal(as.vector(f$query$tra), b$qtra)
    nflags <- nflags + sum(b$tra) + sum(b$qtra)
  }
  expect_gt(nflags, 0)  # The test is not void
})

test_that("a move seen in one genome only is flagged on that side only", {
  # The moved range is alone on its query sequence, so that only the target
  # genome shows a context for it; in the swapped object, only the query
  # genome does.
  f <- flagTranslocations(exampleTranslocation)
  expect_equal(as.vector(f$tra),       c(TRUE,  FALSE, FALSE))
  expect_equal(as.vector(f$query$tra), c(FALSE, FALSE, FALSE))
  f <- flagTranslocations(swap(exampleTranslocation, sort = TRUE))
  expect_equal(as.vector(f$tra),       c(FALSE, FALSE, FALSE))
  expect_equal(as.vector(f$query$tra), c(TRUE,  FALSE, FALSE))
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
