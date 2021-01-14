#' Data Parse
#' @description This function takes in dose-response data, performs data format checks, and 
#' formats the data to a named list, with each list object containing the response data for 
#' a single dose
#' @param dose_response_data A table containing dose-response data from a single study assay. 
#' The data should have two numeric columns with the first containing dose and the second
#' @param dr_threshold A value between 0 and 1. Default 0.8. The percentage of doses that 
#' have at least 3 responses. Allows cases where some doses have few samples. 
#' containing response. 
#' @keywords internal
#' @noRd

# Process dose-response data ----- 
parse_data <- function(dose_response_data, dr_threshold = 0.8) {
  # Data should contain only two columns 
  if (ncol(dose_response_data) != 2) {
    stop("Data should contain 2 columns")
  }
  
  # Doses and responses should be numeric
  if (any(apply(dose_response_data, 2, function(x) !is.numeric(x)))) {
    stop("Doses and responses should be numeric")
  }
  
  # Remove zero-doses
  dose_response_data <- dose_response_data[dose_response_data[,1] != 0,]
  
  # Need at least 4 doses 
  n_doses <- length(unique(dose_response_data[,1]))
  if (n_doses < 4) {
    stop("At least 4 doses required for spline interpolation")
  }
  
  # Want doses to have more than one response
  # Number of doses
  dose_check <- ceiling(n_doses * dr_threshold)
  if (sum(table(dose_response_data[,1]) >= 3) < dose_check) {
    stop(paste0("Insufficient Data. At least ", 
                dose_check, 
                " doses require at least 3 responses"))
  }
  
  # Rename columns 
  colnames(dose_response_data) <- c("dose", "response")
  
  # Reorder rows
  dose_response_data <- dose_response_data[order(dose_response_data$dose),]
  
  # Convert doses to log10 scale
  dose_response_data$log10_dose <- log10(dose_response_data$dose)
  dose_response_data <- dose_response_data[,c("dose", "log10_dose", "response")]
  # Format data as list
  dose_response_data <- split(dose_response_data, dose_response_data[,1])
  dose_response_data
}

