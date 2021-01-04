context("calculate_pod_quantiles.R")

test_that("Are the predicted PODs as expected?", {
  set.seed(471833)
  temp_1 <- calculate_pod_quantiles(test_data_tnt, resample_size = 5)
  expect_length(temp_1, 2)
  expect_length(temp_1[[1]], 5)
  expect_length(temp_1[[2]], 3)
  expect_true(round(temp_1[[2]][2], 4) == 0.4388)
  
  temp_2 <- calculate_pod_quantiles(test_data_toxpod, resample_size = 20)
  expect_length(temp_2[[2]], 3)
  expect_length(temp_2[[1]], 20)
  
  temp_3 <- calculate_pod_quantiles(test_data_toxpod, resample_size = 20, quantile_probs = c(0.25, 0.75))
  expect_length(temp_3[[2]], 2)
  expect_length(temp_3[[1]], 20)

  expect_gt(temp_2[[2]][2], temp_3[[2]][1])
  expect_lt(temp_2[[2]][2], temp_3[[2]][2])
})
