#' Point of Departure Estimation
#' @description
#' Returns quantiles from estimated Point of Departure distribution
#' @param dose_response_data Data frame with rows representing assays/genes and columns representing samples.
#'  The first column should contain assay/gene identifiers. Column names should be numeric values of the
#'  sample's concentration.
#' @param resample_size The number of resamples to perform
#' @param interpolation_size The number of doses to perform spline predictions for response
#' @param quantile_probs The quantiles to report for POD predictions
#' @param single_threshold Not in use yet.
#' @param n_cores The number of cores to use for resampling.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom dplyr bind_rows
#' @importFrom rlang set_names
#' @importFrom rlang .data
#' @importFrom stats predict
#' @importFrom stats quantile
#' @importFrom tidyr unnest_wider
#' @importFrom parallel mclapply
#' @importFrom parallel detectCores
#' @importFrom purrr map
#' @importFrom purrr map_dbl
#' @export

calculate_pod_quantiles <- function(dose_response_data,
                                    resample_size = 1000,
                                    interpolation_size = 50,
                                    quantile_probs = c(0.05, 0.5, 0.95),
                                    single_threshold = 0.5,
                                    n_cores = NULL) {
  # Check number of assays
  dose_response_parsed <- dose_response_data %>% parse_data()
  assay_check <- dose_response_parsed %>%
    group_by(.data$id) %>%
    mutate(n_obs = purrr::map_dbl(.data$data, length)) %>%
    summarize(n_single = mean(.data$n_obs == 1))

  # Remove entire assay if per-concentration single observations surpass given threshold
  single_obs_assays <- assay_check %>%
    filter(.data$n_single >= single_threshold) %>%
    pull(.data$id)

  if (length(single_obs_assays) > 0) {
    message(paste0(
      "Need at least 2 observations for resampling. Removing assays: ",
      paste(single_obs_assays, collapse = ", "),
      "."
    ))
    dose_response_parsed <- dose_response_parsed %>% filter(!.data$id %in% single_obs_assays)
    if (nrow(dose_response_parsed) == 0) {
      stop("Not enough data to perform estimation.")
    }
  }

  # Remove assays with fewer than 4 doses
  # (1) Two points would form a straight line
  # (2) At least four points needed to build spline function
  n_dose_check <- dose_response_parsed %>%
    group_by(.data$id) %>%
    summarize(n_dose = length(.data$dose)) %>%
    filter(.data$n_dose < 4) %>%
    pull(.data$id)

  if (length(n_dose_check) > 0) {
    message(paste0(
      "Need at least 4 doses to build spline model. Removing assays: ",
      paste(n_dose_check, collapse = ", "),
      "."
    ))
    dose_response_parsed <- dose_response_parsed %>% filter(!.data$id %in% n_dose_check)
    if (nrow(dose_response_parsed) == 0) {
      stop("Not enough data to perform estimation.")
    }
  }

  # If the proportion of singletons were fewer than threshold,
  # edit the singletons to be a vector of identical values
  # so that resampling occurs as expected.
  dose_response_parsed <- dose_response_parsed %>%
    mutate(data = map(.data$data, function(vec) if (length(vec) == 1) c(vec, vec) else vec))

  # Set number of cores
  if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1

  # Get unique assays
  assays <- unique(dose_response_parsed$id)
  # pod_quantile_estimates <- setNames(vector("list", length(assays)), assays)
  # pod_quantile_estimates <- lapply(assays, function(i) {
  pod_quantile_estimates <- parallel::mclapply(assays, function(i) {
    dose_response_single <- dose_response_parsed %>% filter(.data$id == i)
    doses <- dose_response_single %>% pull(.data$dose)
    doses_interpolated <- calculate_interpolated_doses(doses)
    pod_resamples <- dose_response_single %>%
      mutate(
        random_responses = map(.data$data, sample, size = resample_size, replace = T),
        random_responses = map(.data$random_responses, set_names, paste0("sample", 1:resample_size))
      ) %>%
      unnest_wider(.data$random_responses) %>%
      select(-.data$id, -.data$dose, -.data$data) %>%
      map(function(x) splines::interpSpline(doses, x)) %>%
      map(function(x) predict(x, doses_interpolated)) %>%
      map(function(pred_obj) calculate_pod_from_menger_curvature(pred_obj)) %>%
      unlist() %>%
      quantile(probs = quantile_probs)
    cat("assay:", i, "done\n")
    data.frame(assay = i, t(pod_resamples), check.names = F)
  }, mc.cores = n_cores)

  suppressWarnings(pod_quantile_estimates %>% bind_rows())
}
