context("MySums")


test_that("mySums function should be consistent", {
  count <- c(1,4,6,0,1)
  expected_count <- c(sum(count[2:5]), sum(count[3:5]), sum(count[4:5]), count[5], 0)
  expect_equal(mySums(count), expected_count)
})