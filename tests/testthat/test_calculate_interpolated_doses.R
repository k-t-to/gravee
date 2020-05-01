context("calculate_interpolated_doses.R")

test_that("Does dose-interpolation behave as expected?", {
  expect_error(calculate_interpolated_doses(1))
  expect_error(calculate_interpolated_doses(c(10, 10, 10)))
  expect_equal(length(calculate_interpolated_doses(c(2, 10))), 50)
  expect_equal(length(calculate_interpolated_doses(c(2, 10), 25)), 25)
})
