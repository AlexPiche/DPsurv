context("Brier score funtion")


test_that("Brier Score function should be consistent", {
  prob <- c(0,0.25,0.5,0.75,1)
  status <- c(0,1,0,1,1)
  expected_result <- mean((prob-status)^2)
  expect_equal(BrierScore(prob, status), expected_result)
})