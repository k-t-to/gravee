context("calculate_pod_quantiles.R")

test_that("Are assays with insufficient data removed?", {
  expect_message(
    calculate_pod_quantiles(test_data_removes, resample_size = 1, n_cores = 1),
    "Need at least 2 observations for resampling."
  )
  expect_message(
    calculate_pod_quantiles(test_data_doses, resample_size = 1, n_cores = 1),
    "Need at least 4 doses to build spline model."
  )
})

test_that("If all assays are removed, does the analysis stop?", {
  expect_error(
    calculate_pod_quantiles(test_data_1, resample_size = 1, n_cores = 1),
    "Not enough data to perform estimation"
  )
  expect_error(
    calculate_pod_quantiles(test_data_doses_error, resample_size = 1, n_cores = 1),
    "Not enough data to perform estimation."
  )
})

test_that("Are the predicted PODs as expected?", {
  temp_1 <- calculate_pod_quantiles(test_data_duplicate, resample_size = 1, n_cores = 1)
  expect_equal(dim(temp_1), c(1, 4))
  expect_equal(as.numeric(round(temp_1[3], 2)), 0.41)

  temp_2 <- calculate_pod_quantiles(test_data_more_duplicate, resample_size = 1, n_cores = 1)
  expect_equal(dim(temp_2), c(2, 4))
  expect_equal(as.numeric(round(temp_2[1, 3], 2)), 0.41)
  expect_equal(as.numeric(round(temp_2[2, 3], 2)), 3.16)

  temp_3 <- calculate_pod_quantiles(test_data_reformat, resample_size = 1, n_cores = 1)
  expect_equal(as.numeric(round(temp_3[1, 3], 2)), 173.02)
  expect_equal(as.numeric(round(temp_3[2, 3], 2)), 211.47)
})
