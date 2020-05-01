#' Menger Curvature
#' @description Given three xy-coordinates, calculates their Menger Curvature.
#' @param interpolated_dose_vector Numeric vector of length 3 containing
#'  doses corresponding to x-axis values.
#' @param predicted_response_vector Numeric vector of length 3 containing
#'  responses corresponding to y-axis values
#' @importFrom stats dist
#' @export

# Calculate Menger Curvature -----
calculate_menger_curvature <- function(interpolated_dose_vector, predicted_response_vector) {
  if (length(interpolated_dose_vector) != 3 | length(predicted_response_vector) != 3) {
    stop("Need 3 data points")
  }

  # Area calculation for numerator
  matrix_temp <- cbind(interpolated_dose_vector, predicted_response_vector, 1)
  A <- 0.5 * abs(det(matrix_temp))
  numerator <- 4 * A

  # Distance calculation for denominator
  line_1 <- dist(matrix_temp[c(1, 2), -3])
  line_2 <- dist(matrix_temp[c(2, 3), -3])
  line_3 <- dist(matrix_temp[c(1, 3), -3])
  # denominator <- prod(
  #  abs(line_1 - line_2),
  #  abs(line_1 - line_3),
  #  abs(line_2 - line_3)
  # )
  denominator <- line_1 * line_2 * line_3

  as.numeric(numerator / denominator)
}
