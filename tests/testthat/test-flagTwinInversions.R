twinRev <- exampleTwinInversions |> reverse() |> sort(ignore.strand = TRUE)

test_that("twin inversions are flagged on their first range", {
  expect_equal(as.vector(flagTwinInversions(exampleTwinInversions)$twi),
               c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(as.vector(flagTwinInversions(twinRev)$twi),
               c(TRUE, FALSE, FALSE, FALSE))
})

test_that("the counter-example is not flagged", {
  expect_false(any(flagTwinInversions(exampleNotTwinInversions)$twi))
})

test_that("ranges on different target sequences are not flagged", {
  for (gb in list(exampleTwinInversions, twinRev)) {
    expect_false(any(flagTwinInversions(moveToOtherTargetSeq(gb, 1))$twi))
    expect_false(any(flagTwinInversions(moveToOtherTargetSeq(gb, 4))$twi))
  }
})

test_that("short objects return a flag of the right length", {
  for (n in 0:3)
    expect_length(flagTwinInversions(exampleTwinInversions[seq_len(n)])$twi, n)
})
