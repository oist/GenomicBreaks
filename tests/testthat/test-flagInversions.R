invRev <- exampleInversion |> reverse() |> sort(ignore.strand = TRUE)

test_that("inversions are flagged on their first range", {
  expect_equal(as.vector(flagInversions(exampleInversion)$inv), c(TRUE, FALSE, FALSE))
  expect_equal(as.vector(flagInversions(invRev)$inv),           c(TRUE, FALSE, FALSE))
})

# Range 1 lies on the query between the first two ranges of a triplet that
# would otherwise be an inversion, so the triplet is not adjacent on the query.
rangeOneInside <- list(
  plus  = GBreaks(target = GRanges(c("chrA:1-50:+", "chrA:100-110:+", "chrA:120-130:-", "chrA:140-150:+")),
                  query  = GRanges(c("chrB:112-118", "chrB:100-110",  "chrB:120-130",   "chrB:140-150"))),
  minus = GBreaks(target = GRanges(c("chrA:1-50:+", "chrA:100-110:-", "chrA:120-130:+", "chrA:140-150:-")),
                  query  = GRanges(c("chrB:132-138", "chrB:140-150",  "chrB:120-130",   "chrB:100-110"))))

test_that("a range between the triplet ranges on the query prevents the flag, even if it is range 1", {
  for (gb in rangeOneInside) {
    expect_false(any(flagInversions(gb)$inv))
    # Same object with an unrelated range in front, so that the range between
    # the triplet ranges is range 2.
    extra <- GBreaks(target = GRanges("chrA:1-5:+"), query = GRanges("chrC:1-5"))
    seqlevels(extra$query) <- c("chrB", "chrC")
    gb2 <- c(extra, gb)
    expect_true(isSorted(gb2))
    expect_false(any(flagInversions(gb2)$inv))
  }
})
