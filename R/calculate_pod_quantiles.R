#' Point of Departure Estimation
#' @description
#' Estimate a chemical Point of Departure (POD) from dose-response data using
#' bootstrap resampling and spline meta-regression
#' @param dose_response_data Data frame with two numeric columns. The first
#' column should contain doses and the second column should contain the responses.  
#' @param resample_size The number of bootstrap samples  
#' @param interpolation_size The number of doses to interpolate between the
#' minimum and maximum doses for spline prediction
#' @param quantile_probs The quantiles to report for POD predictions
#' @return 
#' The returned object prints the POD quantiles to the console. The returned object
#' is a gravee list with four named components. The first component `pods` is a 
#' vector of length `resample_size` containing the estimated PODs from all bootstrap
#' samples, on the original dose scale. The second element `pod_quantiles` contains
#' the POD quantiles, with boundaries defined by `quantile_probs`. The third
#' component `all_res` is a list of length `resample_size` containing the 
#' interpolated spline predictions from each bootstrap sample and curvature
#' measurements along the fit curve. The fourth component `bs_samples` is a list
#' of length `resample_size` containing the sampled values for each bootstrap
#' sample. 
#' @details
#' Analysis is performed on log10-transformed doses. The returned PODs are 
#' back transformed to the original dose scale. 
#' 
#' If the input data contains a dose = 0, the 0 dose is converted to be 1/10th 
#' the minimum non-zero dose, such that after log10-transformation, the distance
#' between the 0 dose and minimum non-zero dose is 1. For PODs and interpolated 
#' doses less than the minimum non-zero dose, values are back-transformed and 
#' rescaled to be left bound by 0. Thus, in the output, for interpolated doses 
#' less than the minimum non-zero dose, values in the `log10_dose` column 
#' will not be equal to log10 of the `dose` column. 
#' 
#' POD estimation assumes that the input dose-response data are asymptotic at 
#' low-doses. If the lower limit of the POD credible interval is within the two
#' smallest doses, the function will return a warning: 
#' `POD Credible Interval violates low-dose asymptote assumption.` In this case,
#' verify that an asymptote has been established at lower doses of the input data. 
#' 
#' @importFrom stats quantile approxfun density median
#' @importFrom graphics axis box grid layout legend par polygon segments
#' @export
#' 
#' @examples
#' x <- rep(c(0,0.05,0.1,0.5,1), each = 3)
#' y <- 1/exp(-x) + abs(rnorm(length(x)))
#' df <- data.frame(x,y)
#' pods <- calculate_pod_quantiles(df)
#' plot(pods)

calculate_pod_quantiles <- function(dose_response_data,
                                    resample_size = 1000,
                                    interpolation_size = 50,
                                    quantile_probs = c(0.025, 0.5, 0.975)) {
  # Parse data to list
  dose_response_parsed <- parse_data(dose_response_data)
  # Get doses
  interpolated_doses <- calculate_interpolated_doses(dose_response_parsed,
                                                     interpolation_size = interpolation_size)
  
  # Perform bootstrap
  mc_res <- lapply(1:resample_size, function(x) {
    perform_bootstrap(dose_response_parsed,
                      interpolated_doses,
                      bs_id = x)
  })

  # Pull PODs
  pods <- lapply(mc_res, function(x) x[[1]][which.max(x[[1]]$mc),])
  pods <- do.call("rbind", pods)
  colnames(pods)[2:3] <- c("pod", "log10_pod")
  rownames(pods) <- NULL

  # Calculate POD quantiles
  pod_quants <- quantile(pods$pod, sort(quantile_probs))

  # Issue warning if Lower Quantile is within the first two doses
  pod_check <- sort(as.numeric(names(dose_response_parsed)))[2]
  if (pod_quants[1] < pod_check) {
    warning("POD Credible Interval violates low-dose asymptote assumption.")
  }
  
  all_mc <- lapply(mc_res, function(x) x[[1]])
  bs_samples <- lapply(mc_res, function(x) x[[2]])

  out <- list(pods = pods$pod,
              pod_quantiles = pod_quants,
              all_res = all_mc,
              bs_samples = bs_samples)
  class(out) <- c("gravee", class(out))
  return(out)
}

# class specific - show only the pod quantiles
#' @export
print.gravee <- function(x, ...) {
  x <- x[[2]]
  print(x, ...)
}

# plot pod distribution
#' @export
plot.gravee <- function(x,
                        xlab = NULL,
                        main = "Distribution of PODs",
                        ...) {
  op_def <- par(no.readonly = T)

  # Get x-values for x-axis scale
  x_vals <- x$bs_samples[[1]][,c("dose", "log10_dose")]
  x_vals$dose <- ifelse(x_vals$dose < 0.01, signif(x_vals$dose, 2), round(x_vals$dose, 2))

  # Get PODs on log scale
  pods <- unlist(lapply(x$all_res, function(x) x$log10_dose[which.max(x$mc)]))
  # Calculate density and density function
  d <- density(pods)
  dd <- approxfun(d)

  # Credible Interval
  # Get legend labels
  quantile_probs <- as.numeric(gsub("%", "", names(x$pod_quantiles)))/100
  quantile_probs <- quantile_probs[c(1, length(quantile_probs))]
  ci <- quantile(x$pods, quantile_probs)
  ci <- ifelse(ci < 0.01, signif(ci, 0.01), round(ci, 2))
  ci_lab <- paste0(quantile_probs*100, "% = ", ci)

  # Get plot line locations
  qls <- quantile(pods, quantile_probs)
  # Get shading region and ymax values
  ql_id <- min(which(d$x >= qls[1]))
  qu_id <- max(which(d$x < qls[2]))
  ql_y1 <- d$y[ql_id]
  qu_y1 <- d$y[qu_id]

  # Median
  # Get legend label
  med <- median(x$pods)
  med <- ifelse(med < 0.01, signif(med, 2), round(med, 2))
  med_lab <- paste0("Median = ", med)
  # Plot line location 
  med_line <- median(pods)

  # Set plot values
  if (is.null(xlab)) xlab <- expression("POD Estimate ("~Log[10]*"Scale)")
  
  layout(matrix(c(1,1,1,1,2), nrow = 1))
  par(cex = 1, mar = c(5.1,4.1,2.1,1.1))
  # Draw base density plot
  plot(d,
       panel.first = grid(),
       xlab = xlab,
       xlim = range(x_vals$log10_dose),
       main = main,
       frame.plot = F,
       xaxt = "n",
       ...)
  axis(1, at = x_vals$log10_dose, labels = x_vals$dose)
  box(lwd = 2, bty = "l")

  # Fill in credible interval
  polygon(d, col = "gray65")
  segments(x0 = qls, x1 = qls, y0 = 0, y1 = c(ql_y1, qu_y1), lwd = 2, col = "blue")
  with(d, polygon(x = x[c(ql_id, ql_id:qu_id, qu_id)],
                  y = c(0, y[ql_id:qu_id], 0),
                  col = "gray25",
                  border = NA))

  # Add median
  segments(x0 = med_line, x1 = med_line, y0 = 0, y1 = dd(med_line),
           lwd = 2, col = "orange")

  # Add legend
  par(cex = 0.65, mar = c(5.1,1.1,2.1,1.1))
  plot(0, 0, axes = F, ann = F, type = "n", xlim = c(-1,1), ylim = c(-1,1))
  legend(-1.25,0.25, legend = ci_lab, col = "blue", lty = "solid", lwd = 2, title = "POD Credible Interval", xpd = T)
  legend(-1.25,-0.15, legend = med_lab, col = "orange", lty = "solid", lwd = 2, xpd = T)
  par(op_def)
  layout(1)
}
