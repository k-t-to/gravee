#' Point of Departure Estimation
#' @description
#' Returns quantiles from estimated Point of Departure distribution
#' @param dose_response_data Data frame with rows representing assays/genes and columns representing samples.
#'  The first column should contain assay/gene identifiers. Column names should be numeric values of the
#'  sample's concentration.
#' @param resample_size The number of resamples to perform
#' @param interpolation_size The number of doses to perform spline predictions for response
#' @param quantile_probs The quantiles to report for POD predictions
#' @param dr_threshold A value between 0 and 1. Default 0.8. The percentage of doses that 
#' have at least 3 responses. Allows cases where some doses have few samples. 
#' @importFrom stats quantile
#' @export

calculate_pod_quantiles <- function(dose_response_data,
                                    resample_size = 1000,
                                    interpolation_size = 50,
                                    quantile_probs = c(0.05, 0.5, 0.95),
                                    dr_threshold = 0.8) {
  # Parse data to list
  dose_response_parsed <- parse_data(dose_response_data, dr_threshold)
  # Get doses
  dose_vector <- sapply(dose_response_parsed, function(x) x$log10_dose[1])
  interpolated_doses <- seq(min(dose_vector), max(dose_vector), length.out = interpolation_size)
  
  # Perform bootstrap
  pods <- lapply(1:resample_size, function(x) perform_bootstrap(dose_response_parsed, interpolated_doses))
  
  # Calculate quantiles
  pods_dose <- sapply(pods, function(x) x["pod"])
  names(pods_dose) <- NULL
  pod_quants <- quantile(pods_dose, quantile_probs)
  # Return pod vector and quantiles
  list(
    pod_distribution = pods_dose,
    pod_quantiles = pod_quants
  )
}
