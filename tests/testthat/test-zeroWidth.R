# Zero-width ranges, such as the mapped regions between adjacent alignments.

test_that("zero-width bridges are colinear with their flanks", {
  for (gb in list(withBridges(exampleColinear3),
                  withStrandedBridges(exampleColinear3),
                  withStrandedBridges(reverse(exampleColinear3, query = TRUE) |> sort(ignore.strand = TRUE)))) {
    expect_true(any(width(gb) == 0))
    expect_equal(as.vector(flagColinearAlignments(gb)$colinear), c(TRUE, TRUE, TRUE, TRUE, FALSE))
    expect_equal(as.character(flagPairs(gb)$pairs), c("Col", "Col", "Col", "Col", "Bnd"))
    f <- flagTranslocations(gb)
    expect_false(any(f$tra) || any(f$trs))
    expect_length(coalesce_contigs(gb), 1)
  }
})

test_that("mapped regions of width zero on one genome are colinear with their flanks", {
  zeroOnQuery  <- GBreaks(target = GRanges(c("chrA:100-189:+", "chrA:190-199:*", "chrA:200-289:+")),
                          query  = GRanges(c("chrB:100-189",   "chrB:190-189",   "chrB:190-279")))
  zeroOnTarget <- GBreaks(target = GRanges(c("chrA:100-189:+", "chrA:190-189:*", "chrA:190-279:+")),
                          query  = GRanges(c("chrB:100-189",   "chrB:190-199",   "chrB:200-289")))
  for (gb in list(zeroOnQuery, zeroOnTarget)) {
    expect_equal(as.vector(flagColinearAlignments(gb)$colinear), c(TRUE, TRUE, FALSE))
    expect_equal(as.character(flagPairs(gb)$pairs), c("Col", "Col", "Bnd"))
  }
})

test_that("zero-width bridges do not change the flags of the alignments", {
  flags <- function(gb) {
    f <- flagTranslocations(gb)
    data.frame(col = as.vector(flagColinearAlignments(gb)$colinear),
               inv = as.vector(flagInversions(gb)$inv),
               dbl = as.vector(flagDoubleInversions(gb)$Dbl),
               twi = as.vector(flagTwinInversions(gb)$twi),
               tra = as.vector(f$tra),
               trs = as.vector(f$trs),
               pairs = as.character(flagPairs(gb)$pairs))
  }
  for (gb in loadSimulations(8)) {
    gb <- adjacent(gb)
    ref <- flags(gb)
    for (b in list(withBridges(gb), withStrandedBridges(gb))) {
      isBridge <- width(b) == 0
      expect_equal(sum(!isBridge), length(gb))
      fb <- flags(b)
      expect_equal(fb[!isBridge, ], ref, ignore_attr = TRUE)
      # Bridges are colinear with the next range, and carry no other flag.
      expect_true(all(fb$col[isBridge]))
      expect_true(all(fb$pairs[isBridge] == "Col"))
      expect_false(any(as.matrix(fb[isBridge, c("inv", "dbl", "twi", "tra", "trs")])))
      expect_equal(as.character(coalesce_contigs(b)), as.character(coalesce_contigs(gb)))
    }
  }
})

test_that("precedeZW() and followZW() equal precede() and follow() without zero-width ranges", {
  for (gb in loadSimulations(10)) {
    expect_identical(precedeZW(gb), precede(gb))
    expect_identical(followZW(gb),  follow(gb))
    expect_identical(precedeZW(gb$query, ignore.strand = TRUE), precede(gb$query, ignore.strand = TRUE))
  }
})
