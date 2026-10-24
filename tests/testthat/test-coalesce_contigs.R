same <- function(a, b) {
  expect_equal(as.character(granges(a)), as.character(granges(b)))
  expect_equal(as.character(a$query),    as.character(b$query))
}

test_that("coalesce_contigs() coalesces over mapped regions", {
  # Bridge regions of width 9
  same(coalesce_contigs(withBridges(exampleColinear5)), coalesce_contigs(exampleColinear5))
  # Bridge regions of width zero on both genomes
  same(coalesce_contigs(withBridges(exampleColinear3)), coalesce_contigs(exampleColinear3))
  # Two colinear runs per strand around an inversion
  gb <- permutationVectorToGBreaks(c(1, 2, -4, -3, 5, 6))
  same(coalesce_contigs(withBridges(gb)), coalesce_contigs(gb))
  expect_length(coalesce_contigs(gb), 3)
})

test_that("a mapped region of width zero on one genome is absorbed", {
  # The target genome has an unaligned region (for instance a new intron)
  # between two alignments that are adjacent on the query genome.
  gb <- GBreaks(target = GRanges(c("chrA:100-189:+", "chrA:190-199:*", "chrA:200-289:+")),
                query  = GRanges(c("chrB:100-189",   "chrB:190-189",   "chrB:190-279")))
  coa <- coalesce_contigs(gb)
  expect_equal(as.character(granges(coa)), "chrA:100-289:+")
  expect_equal(as.character(coa$query),    "chrB:100-279")
})

test_that("mapped regions outside of coalesced ranges are kept", {
  gb <- exampleColinear3
  extra <- GBreaks(target = GRanges("chrA:450-460:*"), query = GRanges("chrB:500-510"))
  seqinfo(extra) <- seqinfo(gb)
  seqinfo(extra$query) <- seqinfo(gb$query)
  coa <- coalesce_contigs(c(gb, extra))
  expect_equal(as.character(granges(coa)), c("chrA:100-400:+", "chrA:450-460:*"))
  expect_equal(score(coa), width(coa))
})

test_that("objects where all ranges are unstranded are coalesced", {
  expect_length(coalesce_contigs(unstrand(exampleColinear3)), 1)
  expect_length(coalesce_contigs(unstrand(exampleInversion)), 1)
})
