#' General functions
#' @description Contains general helper functions
#' @param dose_vector numeric vector of doses on the log scale
#' @param min_dose minimum non zero dose (original scale)
#' @keywords internal
#' @noRd

# Rescale values between two numbers -----
rescale <- function(original, old_min = NULL, old_max = NULL, new_min, new_max) {
  if (is.null(old_min)) old_min <- min(original)
  if (is.null(old_max)) old_max <- max(original)
  
  # Append old min and max to ensure proper scaling
  original <- c(old_min, old_max, original)
  
  scale_val <- (new_max - new_min)/(old_max - old_min)
  new_values <- (scale_val * (original - old_min)) + new_min
  
  # Remove appended values
  new_values[-c(1:2)]
}

# Backtransform small values -----
backtransform_small <- function(dose_vector, min_dose) {
  dose_transform <- 10^dose_vector
  # Which values are smaller than the minimum non-zero dose?
  is_small <- dose_transform <= min_dose
  if (any(is_small)) {
    dose_small <- dose_transform[is_small]
    dose_small <- rescale(dose_small,
                          old_min = min_dose/10,
                          old_max = min_dose,
                          new_min = 0,
                          new_max = min_dose)
    dose_transform[is_small] <- dose_small
  }
  dose_transform
}

# Calculate interpolated doses -----
calculate_interpolated_doses <- function(dose_response_parsed,
                                         interpolation_size = 50) {
  # Get log transformed doses
  dose_vector <- sapply(dose_response_parsed, function(x) x$log10_dose[1])
  interpolated_doses_log <- seq(min(dose_vector),
                            max(dose_vector),
                            length.out = interpolation_size)
  # Back transform interpolated doses
  doses <- sort(as.numeric(names(dose_response_parsed)))
  min_dose <- ifelse(doses[1] == 0, doses[2], 0)
  dose_transform <- backtransform_small(interpolated_doses_log, min_dose)
  data.frame(dose = dose_transform,
             log10_dose = interpolated_doses_log)
}
