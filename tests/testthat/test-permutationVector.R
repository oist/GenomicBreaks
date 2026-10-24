test_that("permutationVector() needs oriented ranges", {
  expect_equal(permutationVector(exampleInversion), c(1L, -2L, 3L))
  expect_error(permutationVector(unstrand(exampleInversion)), "orientation")
  expect_error(permutationVector(withBridges(exampleColinear5)), "orientation")
})

test_that("inversionDistance() and breakpointGraphProperties() stop on unstranded ranges", {
  expect_equal(inversionDistance(exampleInversion), 1)
  expect_error(inversionDistance(unstrand(exampleInversion)), "orientation")
  expect_error(breakpointGraphProperties(unstrand(exampleInversion)), "orientation")
})
