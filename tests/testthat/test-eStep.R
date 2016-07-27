context("eStep funtion")

DataStorage <- new("DataStorage")
DataStorage@computation <- c(0, 0)
DataStorage@censoring <- c(0, 1)

test_that("eStep probability must be accurate for censored and uncensored data", {
  log_prob <- eStep(0, 1, 1, DataStorage)
  expect_equal(pnorm(DataStorage@computation[1]), exp(log_prob[1]))
  expect_equal(dnorm(DataStorage@computation[2]), exp(log_prob[2]))
})