#' Data Parse
#' @description This function takes in dose-response data, performs data format checks, and 
#' formats the data to a named list, with each list object containing the response data for 
#' a single dose
#' @param dose_response_data Data frame with two numeric columns. The first column should contain 
#' doses and the second column should contain the responses.  
#' @keywords internal
#' @noRd

# Process dose-response data ----- 
parse_data <- function(dose_response_data) {
  # Data should contain only two columns
  if (ncol(dose_response_data) != 2) {
    stop("Data should contain 2 columns")
  }

  # Remove rows with missing data
  dose_response_data <- dose_response_data[!is.na(dose_response_data[,1]),]
  dose_response_data <- dose_response_data[!is.na(dose_response_data[,2]),]

  # Doses and responses should be numeric
  if (any(apply(dose_response_data, 2, function(x) !is.numeric(x)))) {
    stop("Doses and responses should be numeric")
  }

  # Remove doses with fewer than three replicates
  dose_rm <- as.numeric(names(which(table(dose_response_data[,1]) < 3)))
  if (length(dose_rm) > 0) {
    warning(paste0("Minimum of 3 replicates per dose. Removing doses from dataset: ",
                   paste(dose_rm, collapse = ", ")))
    dose_response_data <- dose_response_data[!dose_response_data[,1] %in% dose_rm,]
  }

  # Need at least 4 doses
  n_doses <- length(unique(dose_response_data[,1]))
  if (n_doses < 4) {
    stop("At least 4 doses with 3 replicates required for spline interpolation")
  }

  # Rename columns
  colnames(dose_response_data) <- c("dose", "response")

  # Reorder rows
  dose_response_data <- dose_response_data[order(dose_response_data$dose),]

  # Convert doses to log10 scale
  dose_response_data$log10_dose <- log10(dose_response_data$dose)
  if (any(dose_response_data$dose == 0)) {
    min_val <- min(dose_response_data$dose[dose_response_data$dose != 0])
    new_zero <- log10(min_val/10)
    dose_response_data$log10_dose[dose_response_data$dose == 0] <- new_zero
  } 

  dose_response_data <- dose_response_data[,c("dose", "log10_dose", "response")]
  # Format data as list
  dose_response_data <- split(dose_response_data, dose_response_data[,1])
  dose_response_data
}