test_that("dist2next() returns the number of bases between ranges", {
  gr <- GRanges(c("chrA:100-150", "chrA:151-200", "chrA:190-250", "chrA:261-300", "chrB:1-10"))
  # adjacent, overlapping, 10-base gap, other sequence, last range
  expect_equal(dist2next(gr)$dist, c(0, 0, 10, NA, Inf))
  # 151-189 and 201-260 are between the ranges two steps apart
  expect_equal(dist2next(gr, step = 2)$dist, c(39, 60, NA, Inf, Inf))
})

test_that("dist2next() on GBreaks objects", {
  # exampleInversion: gaps of 19 and 9 bases on both genomes
  d <- dist2next(exampleInversion, ignore.strand = TRUE)
  expect_equal(d$tdist, c(19, 9, Inf))
  expect_equal(d$qdist, c(19, 9, Inf))
})

test_that("dist2next() accepts objects with fewer ranges than step", {
  expect_length(dist2next(exampleInversion[0])$tdist, 0)
  expect_length(dist2next(granges(exampleInversion)[0])$dist, 0)
  expect_equal(dist2next(exampleInversion[1], step = 2)$qdist, Inf)
  expect_equal(dist2next(exampleInversion, step = 4)$tdist, rep(Inf, 3))
})

test_that("dist2next() leaves the metadata of the query ranges unchanged", {
  gb <- exampleInversion
  gb$query$dist <- 1:3
  expect_equal(dist2next(gb)$query$dist, 1:3)
})

test_that("tol bridges gaps of up to tol bases", {
  # exampleColinear: ranges separated by 100 bases on both genomes
  expect_true (flagColinearAlignments(exampleColinear, tol = 100)$colinear[1])
  expect_false(flagColinearAlignments(exampleColinear, tol =  99)$colinear[1])
  expect_length(coalesce_contigs(exampleColinear, tol = 100), 1)
  expect_length(coalesce_contigs(exampleColinear, tol =  99), 2)
  expect_equal(as.vector(flagInversions(exampleInversion, tol = 19)$inv), c(TRUE,  FALSE, FALSE))
  expect_equal(as.vector(flagInversions(exampleInversion, tol = 18)$inv), c(FALSE, FALSE, FALSE))
})

test_that("tol = 0 coalesces adjacent ranges only", {
  adj <- GBreaks(target = GRanges(c("chrA:100-150:+", "chrA:151-200:+", "chrA:202-250:+")),
                 query  = GRanges(c("chrB:100-150",   "chrB:151-200",   "chrB:202-250")))
  coa <- coalesce_contigs(adj, tol = 0)
  expect_equal(as.character(granges(coa)), c("chrA:100-200:+", "chrA:202-250:+"))
  expect_equal(as.character(coa$query),    c("chrB:100-200",   "chrB:202-250"))
})
