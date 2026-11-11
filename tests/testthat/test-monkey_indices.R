yeast <- load_genomic_breaks(
  system.file("extdata/SacCer3__SacPar.gff3.gz", package = "GenomicBreaks"))
yeast$aLength  <- pmax(width(yeast), width(yeast$query))
yeast$gapOpens <- as.integer(width(yeast) != width(yeast$query))

geometric <- list(wukong_index, abu_index, albert_index, marcel_index, imo_index,
                  rafiki_index, hanuman_index, jack_index)

test_that("monkey indices give the reference values", {
  # Values checked against the Python implementation used for the screen.
  expect_equal(wukong_index(yeast),   2.49607603555432)
  expect_equal(abu_index(yeast),      0.122758247423296)
  expect_equal(albert_index(yeast),  -1.08211785824754)
  expect_equal(marcel_index(yeast),   0.128113761372173)
  expect_equal(imo_index(yeast),      0.91634832504211)
  expect_equal(rafiki_index(yeast),   0.039850511875284)
  expect_equal(hanuman_index(yeast),  0.57776027669833)
  expect_equal(jack_index(yeast),    -0.516659528016429)
  expect_equal(george_index(yeast),   0.117812275194519)
  expect_equal(baker_index(yeast),   -0.0141499900111354)
  expect_equal(albert_index(exampleGenomicBreaks), 13.2012180062151)
  expect_equal(jack_index(exampleGenomicBreaks),    0.0684304861068656)
})

test_that("monkey indices do not depend on the representation", {
  y <- forceSeqLengths(yeast)
  for (f in geometric) {
    expect_equal(f(reverse(y)),               f(y))
    expect_equal(f(reverse(y, query = TRUE)), f(y))
    expect_equal(f(y[sample(length(y))]),     f(y))  # order of the ranges
  }
  expect_equal(george_index(reverse(y)), george_index(y))
  expect_equal(baker_index(reverse(y)),  baker_index(y))
  # order of the sequence levels
  z <- y
  seqlevels(z) <- rev(seqlevels(z))
  seqlevels(z$query) <- rev(seqlevels(z$query))
  for (f in geometric) expect_equal(f(z), f(y))
})

test_that("monkey indices handle objects that they can not summarise", {
  for (f in geometric) {
    expect_true(is.na(f(exampleInversion[1:2])))
    expect_error(f(unstrand(exampleInversion)), "orientation")
  }
  expect_error(george_index(exampleGenomicBreaks), "score")
  expect_error(baker_index(exampleGenomicBreaks),  "gapOpens")
})
