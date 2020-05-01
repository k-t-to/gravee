context("calculate_pod_from_menger_curvature.R")

test_that("", {
  sim_predictedDR <- list(
    x = 1:10,
    y = sin(1:10)
  )
  # Does the function return the correct pod estimate?
  expect_equal(calculate_pod_from_menger_curvature(sim_predictedDR), 8)
})
