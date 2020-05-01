#' Interpolate Doses
#' @description Returns sequence of interpolated doses.
#' @param dose_vector Vector of experimental doses.
#' @param length.out Desired length of sequence. Sent to `seq`
#' @noRd

calculate_interpolated_doses <- function(dose_vector, length.out = 50) {
  if (length(dose_vector) < 2) stop("Cannot interpolate with one dose")
  start_value <- min(dose_vector)
  end_value <- max(dose_vector)
  if (start_value == end_value) stop("No min or max found")
  # new_doses <- seq(from = start_value, to = end_value, length.out = length.out + 1)
  new_doses <- seq(from = start_value, to = end_value, length.out = length.out)
  new_doses
  # new_doses[-1]
}
