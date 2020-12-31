context("calculate_pod_quantiles.R")

test_that("Are the predicted PODs as expected?", {
  temp_1 <- calculate_pod_quantiles(test_data_tnt, resample_size = 1)
  expect_length(temp_1, 3)
  expect_true(round(temp_1[2], 4) == 0.4415)
  
  temp_2 <- calculate_pod_quantiles(test_data_toxpod, resample_size = 20)
  temp_3 <- calculate_pod_quantiles(test_data_toxpod, resample_size = 20, quantile_probs = c(0.25, 0.75))
  expect_length(temp_3, 2)
  expect_gt(temp_2[2], temp_3[1])
  expect_lt(temp_2[2], temp_3[2])
})
