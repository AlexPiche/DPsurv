context("Log score funtion")


test_that("Log Score function should be consistent", {
  prob <- c(0,0.25,0.5,0.75,1)
  status <- c(0,1,0,1,1)
  expected_result <- mean(ifelse((1-status),log(1-prob),log(prob)))
  expect_equal(MeanLogScore(prob, status), expected_result)
})