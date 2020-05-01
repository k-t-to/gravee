context("parse_data.R")

test_that("Are the data being prepared correctly?", {
  # Does the function return the correct dimensions?
  expect_equal(dim(parse_data(test_data_logcpm)), c(60, 3))
  # Does the function return the correct data when there is one observation per assay?
  expect_length(unlist(parse_data(test_data_1)$data), 8)
})
