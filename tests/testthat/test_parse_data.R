context("parse_data.R")

test_that("Are the data being prepared correctly?", {
  # Does the function return a list with the correct dimensions
  parse_test_1 <- parse_data(test_data_tnt)
  expect_length(parse_test_1, 6)
  expect_equal(names(parse_test_1), c("0", "0.3125", "0.625", "1.25", "2.5", "5"))
  expect_true(all(matrix(c(4,3), ncol = 6, nrow = 2) == sapply(parse_test_1, dim)))
  
  # Are errors produced? 
  # Number of columns
  expect_error(parse_data(data.frame(test_data_tnt[,1])),
               "Data should contain 2 columns")
  expect_error(parse_data(data.frame(test_data_tnt, test_data_tnt)),
               "Data should contain 2 columns")
  # Non numeric
  expect_error(parse_data(data.frame(a = c("dog", "cat"), b = 1:2)),
               "Doses and responses should be numeric")
  # Number of doses
  expect_error(parse_data(test_data_tnt[test_data_tnt$dose %in% c(1.25, 5),]),
               "At least 4 doses with 3 replicates required for spline interpolation")
  # Doses with insufficient data
  parse_test_2 <- rbind.data.frame(test_data_tnt[!test_data_tnt$dose %in% c(1.25,2.5,5),],
                                   test_data_tnt[test_data_tnt$dose == 5,][1,],
                                   test_data_tnt[test_data_tnt$dose == 2.5,][1,],
                                   test_data_tnt[test_data_tnt$dose == 1.25,][1,])
  expect_error(expect_warning(parse_data(parse_test_2), 
                              "Minimum of 3 replicates per dose. Removing doses from dataset: 1.25, 2.5, 5"),
               "At least 4 doses with 3 replicates required for spline interpolation")
})
