#' Point of Departure Estimation
#' @description
#' Estimate a chemical Point of Departure (POD) from dose-response data using bootstrap resampling
#' and spline meta-regression  
#' @param dose_response_data Data frame with two numeric columns. The first column should contain 
#' doses and the second column should contain the responses.  
#' @param resample_size The number of bootstrap samples  
#' @param interpolation_size The number of doses to interpolate between the minimum and maximum doses for
#' spline prediction
#' @param quantile_probs The quantiles to report for POD predictions
#' @param spline_res Logical, default = False. Return the output of the spline fit
#' @return 
#' A list is returned. The first component is a numeric vector of length `resample_size` containing the 
#' estimated PODs from all bootstrap samples, on the original dose scale. The second element returns
#' the credible interval of PODs with the boundaries defined by `quantile_probs`. If `spline_res=T`, the 
#' third element of the list is a data frame with spline fit predictions for all bootstrap samples. 
#' @details
#' Analysis is performed on log10-transformed doses. The returned PODs and interpolated doses are 
#' back transformed to the original dose scale. 
#' 
#' If the input data contains a dose = 0, the 0 dose is converted to be 1/10th the minimum non-zero
#' dose, such that after log10-transformation, the distance between the 0 dose and minimum non-zero 
#' dose is 1. For PODs and interpolated doses less than the minimum non-zero dose, values are 
#' back-transformed and rescaled to be left bound by 0. Thus, in the spline result data frame, 
#' for interpolated doses less than the minimum non-zero dose, values in the `log10dose` column 
#' will not be equal to log10 of the `dose` column. 
#' 
#' POD estimation assumes that the input dose-response data are asymptotic at low-doses. If the
#' lower limit of the POD credible interval is within the two smallest doses, the function will
#' return a warning: `POD Credible Interval violates low-dose asymptote assumption.` In this 
#' case, verify that an asymptote has been established at lower doses of the input data. 
#' 
#' @importFrom stats quantile
#' @export

calculate_pod_quantiles <- function(dose_response_data,
                                    resample_size = 1000,
                                    interpolation_size = 50,
                                    quantile_probs = c(0.025, 0.5, 0.975),
                                    spline_res = F) {
  # Parse data to list
  dose_response_parsed <- parse_data(dose_response_data)
  # Get doses
  dose_vector <- sapply(dose_response_parsed, function(x) x$log10_dose[1])
  interpolated_doses <- seq(min(dose_vector), max(dose_vector), length.out = interpolation_size)
  
  # Perform bootstrap
  pods <- lapply(1:resample_size, function(x) perform_bootstrap(dose_response_parsed, interpolated_doses, spline_res=spline_res))
  
  # Rescale PODs
  doses <- as.numeric(names(dose_response_parsed))
  pods_dose <- sapply(pods, function(x) x[[1]]["pod"])
  names(pods_dose) <- NULL
  
  if (any(doses == 0)) {
    min_dose <- min(doses[doses != 0])
    is_small <- pods_dose <= min_dose
    if (any(is_small)) {
      pod_small <- pods_dose[is_small]
      trans_0 <- min_dose/10
      pod_small <- c(trans_0, min_dose, pod_small)
      pod_small <- min_dose * ((pod_small - min(pod_small))/(max(pod_small) - min(pod_small)))
      pod_small <- pod_small[-c(1:2)]
      pods_dose[is_small] <- pod_small
    }
  } 
  
  pod_quants <- quantile(pods_dose, quantile_probs)
  
  # Issue warning if Lower Quantile is within the first two doses
  pod_check <- sort(as.numeric(names(dose_response_parsed)))[2]
  if (pod_quants[1] < pod_check) {
    warning("POD Credible Interval violates low-dose asymptote assumption.")
  }
  
  # Create spline output and rescale doses
  if (spline_res) {
    spline_output <- lapply(1:resample_size, function(i){
      data.frame(bs_id = i,
                 log10dose = pods[[i]][[2]]$x,
                 pred_response = pods[[i]][[2]]$y)
    })
    spline_output <- do.call("rbind", spline_output)
    spline_doses <- 10^spline_output$log10dose
    if (any(doses == 0)) {
      is_small <- spline_doses <= min_dose
      trans_0 <- min_dose/10
      spline_small <- c(trans_0, min_dose, spline_doses[is_small])
      spline_small <- min_dose * ((spline_small - min(spline_small))/(max(spline_small) - min(spline_small)))
      spline_small <- spline_small[-c(1:2)]
      spline_doses[is_small] <- spline_small
    }
    spline_output$dose <- spline_doses
    spline_output <- spline_output[,c("bs_id", "dose", "log10dose", "pred_response")]
    out <- list(
      pod_distribution = pods_dose,
      pod_quantiles = pod_quants,
      spline_output = spline_output
    )
  } else {
    out <- list(
      pod_distribution = pods_dose,
      pod_quantiles = pod_quants
    )
  }
  return(out)
}
