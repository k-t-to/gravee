#' Calculate POD from Menger Curvature
#' @description Calculates Menger Curvature across interpolated dose-response
#'  and returns the dose with the highest Menger Curvature.
#' @param predicted_dose_response An object from `predict.npolyspline` of length 2 with
#'  the first object a vector of interpolated doses and the second object a fector
#'  of predicted responses from an interpolation spline.
#' @noRd

# Calculate sequential Menger Curvatures and return the max -----
calculate_pod_from_menger_curvature <- function(predicted_dose_response) {
  # Declare list to store calculations
  n_3 <- length(predicted_dose_response$x) - 2
  MC_values <- list(
    dose = vector("double", n_3),
    curvature = vector("double", n_3)
  )
  # Loop through data and calculate MC
  for (i in 1:n_3) {
    end <- i + 2
    dose_temp <- predicted_dose_response$x[i:end]
    response_temp <- predicted_dose_response$y[i:end]
    MC_temp <- calculate_menger_curvature(
      interpolated_dose_vector = dose_temp,
      predicted_response_vector = response_temp
    )
    MC_values[["dose"]][i] <- dose_temp[2]
    MC_values[["curvature"]][i] <- MC_temp
  }
  # Return the dose corresponding to the highest curvature.
  MC_values$dose[which.max(MC_values$curvature)]
}
