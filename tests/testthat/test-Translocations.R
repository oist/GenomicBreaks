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
