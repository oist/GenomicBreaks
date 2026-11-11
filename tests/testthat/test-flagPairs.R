test_that("pairs of the man page examples", {
  expect_equal(as.character(flagPairs(exampleInversion)$pairs),     c("Inv", "Inv", "Bnd"))
  expect_equal(as.character(flagPairs(exampleTranslocation)$pairs), c("Scr", "Scr", "Bnd"))
})

test_that("range 1 is seen when testing query adjacency", {
  # Range 1 lies on the query between ranges 2 and 3.
  gb <- GBreaks(target = GRanges(c("chrA:1-50:+", "chrA:100-110:+", "chrA:120-130:-", "chrA:140-150:+")),
                query  = GRanges(c("chrB:112-118", "chrB:100-110",  "chrB:120-130",   "chrB:140-150")))
  expect_equal(as.character(flagPairs(gb)$pairs), c("Flp", "Scr", "Inv", "Bnd"))
  # Ranges 4 and 5 of the Garg et al. (2019) example are separated by range 1
  # on the query.
  expect_equal(as.character(flagPairs(exampleInversionGarg2019)$pairs)[4], "Scr")
})
