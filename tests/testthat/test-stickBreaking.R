context("Stick breaking funtion")


test_that("Stick breaking function should be consistent", {
  beta_rv <- c(1/4, 1/3, 5/6, 1/7)
  stick <- c(1/4, 1/3*(1-1/4), 5/6*(1-1/3)*(1-1/4), 1*(1-5/6)*(1-1/3)*(1-1/4))
  expect_equal(stickBreaking(beta_rv), stick)
})