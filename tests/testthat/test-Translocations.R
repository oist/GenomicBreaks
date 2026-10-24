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

test_that("tol limits the distances between adjacent ranges", {
  # exampleTranslocation: U and V are 100 bases apart on the query genome;
  # in the swapped object, X and Y are 100 bases apart on the target genome.
  expect_true (flagTranslocations(exampleTranslocation, tol = 100)$tra[1])
  expect_false(flagTranslocations(exampleTranslocation, tol =  99)$tra[1])
  sw <- swap(exampleTranslocation, sort = TRUE)
  expect_true (flagTranslocations(sw, tol = 100)$trs[1])
  expect_false(flagTranslocations(sw, tol =  99)$trs[1])

  # Target X Y U A V, query X A Y U V, all gaps of 10 bases.  Narrowing A on
  # one genome widens one of its gaps to 60 bases, which only affects the
  # flag for which A must be adjacent to its neighbours on that genome.
  gb <- permutationVectorToGBreaks(c(1, 3, 4, 2, 5))
  narrowQ <- gb; narrowQ$query[4] <- resize(narrowQ$query[4], 40)
  narrowT <- gb; ranges(narrowT)[4] <- resize(ranges(narrowT)[4], 40)
  flags <- function(x, tol) {
    f <- flagTranslocations(x, tol = tol)
    c(tra = as.vector(f$tra)[3], trs = as.vector(f$trs)[1])
  }
  expect_equal(flags(narrowQ, 60), c(tra = TRUE,  trs = TRUE))
  expect_equal(flags(narrowQ, 59), c(tra = TRUE,  trs = FALSE))
  expect_equal(flags(narrowT, 60), c(tra = TRUE,  trs = TRUE))
  expect_equal(flags(narrowT, 59), c(tra = FALSE, trs = TRUE))

  # In the simulations, all gaps between consecutive blocks are 10 bases.
  for (gb in loadSimulations(20)) {
    expect_equal(flagTranslocations(gb, tol = 10), flagTranslocations(gb))
    f <- flagTranslocations(gb, tol = 9)
    expect_false(any(f$tra) || any(f$trs))
  }
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
