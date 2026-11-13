test_that("slidingWindow cuts target and query ranges at window boundaries", {
  # Windows 1-250, 126-375 and 251-500 on the 600-bp target sequence.
  sw <- slidingWindow(exampleInversion, windowSize = 250, stepSize = 125)
  expect_length(sw, 3)
  expect_equal(as.character(sw[[1]]),        c("chrA:100-190:+", "chrA:210-250:-"))
  expect_equal(as.character(sw[[2]]),        c("chrA:126-190:+", "chrA:210-291:-", "chrA:301-375:+"))
  expect_equal(as.character(sw[[3]]),        c("chrA:251-291:-", "chrA:301-400:+"))
  # On the minus strand, the start of the target range matches the end of the
  # query range.
  expect_equal(as.character(sw[[1]]$query),  c("chrB:100-190", "chrB:251-291"))
  expect_equal(as.character(sw[[2]]$query),  c("chrB:126-190", "chrB:210-291", "chrB:301-375"))
  expect_equal(as.character(sw[[3]]$query),  c("chrB:210-250", "chrB:301-400"))
})

test_that("slidingWindow scales scores and can return uncut or merged ranges", {
  gb <- exampleInversion
  gb$score <- c(91, 82, 100)
  sw <- slidingWindow(gb, windowSize = 250, stepSize = 125)
  expect_equal(sw[[1]]$score, c(91, 41))
  uncut <- slidingWindow(gb, windowSize = 250, stepSize = 125, cut = FALSE)
  expect_equal(granges(unlist(uncut, use.names = FALSE)), granges(gb)[c(1, 2, 1, 2, 3, 2, 3)])
  merged <- slidingWindow(gb, windowSize = 250, stepSize = 125, merged = TRUE)
  expect_equal(merged$window_id, c(1, 2, 1, 2, 3, 2, 3))
})
