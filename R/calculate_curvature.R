#' Menger Curvature Calculations
#' @description Functions for calculating Menger Curvature along bootstrapped
#' dose-response curves
#' @param interpolated_dose_vector Numeric vector of length 3 containing
#'  doses corresponding to x-axis values.
#' @param predicted_response_vector Numeric vector of length 3 containing
#'  responses corresponding to y-axis values
#' @param predicted_dose_response An object from `predict.npolyspline` of length 2 with
#'  the first object a vector of interpolated doses and the second object a vector
#'  of predicted responses from an interpolation spline.
#' @param dose_response_parsed input data for main analysis
#' @param interpolated_doses sequence of doses to predict responses from spline model
#' @importFrom stats dist
#' @importFrom stats predict
#' @keywords internal
#' @noRd

# Calculate Menger Curvature -----
# Given three xy-coordinates, calculates their Menger Curvature.
calculate_menger_curvature <- function(interpolated_dose_vector,
                                       predicted_response_vector) {
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

  denominator <- line_1 * line_2 * line_3

  as.numeric(numerator / denominator)
}

# Calculate Menger Curvature given spline predictions
calculate_spline_curvature <- function(predicted_dose_response) {
  # Declare list to store calculations
  n_3 <- length(predicted_dose_response$x) - 2
  MC_values <- list(
    log10_dose = vector("double", n_3),
    mc         = vector("double", n_3)
  )

  # Loop through data and calculate MC
  for (i in 1:n_3) {
    end <- i + 2
    dose_temp <- predicted_dose_response$x[i:end]
    response_temp <- predicted_dose_response$y[i:end]
    MC_temp <- calculate_menger_curvature(
      interpolated_dose_vector  = dose_temp,
      predicted_response_vector = response_temp
    )
    MC_values[["log10_dose"]][i] <- dose_temp[2]
    MC_values[["mc"]][i] <- MC_temp
  }

  # Return curvature for interpolated data
  spline_mc <- data.frame(log10_dose    = predicted_dose_response$x,
                          response_pred = predicted_dose_response$y)
  spline_mc$mc <- MC_values$mc[match(spline_mc$log10_dose, MC_values$log10_dose)]
  spline_mc
}

# Perform bootstrap sampling and curve calculations ----- 
perform_bootstrap <- function(dose_response_parsed, interpolated_doses, bs_id) {
  bootstrap_responses <- lapply(dose_response_parsed,
                                function(x) x[sample(nrow(x),1),])
  bootstrap_responses <- do.call("rbind", bootstrap_responses)

  # Create spline model
  bs_spline_model <- splines::interpSpline(bootstrap_responses[,2],
                                           bootstrap_responses[,3])

  # Predict responses for interpolated doses
  pred_vals <- predict(bs_spline_model, interpolated_doses$log10_dose)

  # Calculate Menger Curvature
  mc <- calculate_spline_curvature(pred_vals)
  # Add back original doses
  mc$dose <- interpolated_doses$dose[match(mc$log10_dose, interpolated_doses$log10_dose)]
  # Add bootstrap id
  mc$bs_index <- bs_id
  mc <- mc[,c("bs_index", "dose", "log10_dose", "response_pred", "mc")]
  bootstrap_responses$bs_index <- bs_id
  bootstrap_responses <- bootstrap_responses[,c("bs_index", "dose", "log10_dose", "response")]

  # Return menger curvature dataframe and bootstrap sample dataframe
  list(mc = mc,
       bootstrap_sample = bootstrap_responses)
}

