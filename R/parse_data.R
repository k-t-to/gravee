#' Data Parse
#' @description This function takes in dose-response data and processes it to a nested tibble
#' @param dose_response_data A table with rows representing assays/genes and columns representing
#'  samples. The first column should contain assay/gene identifiers. Column names should be numeric
#'  values of the sample's concentration.
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom dplyr ungroup
#' @importFrom dplyr mutate
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @noRd

# Parse dose response data to nested tibble -----
parse_data <- function(dose_response_data) {
  names(dose_response_data)[1] <- "id"
  # Create nested dataset
  suppressWarnings(
    dose_response_data %>%
      pivot_longer(cols = -1, names_to = "dose") %>%
      filter(!is.na(.data$value)) %>%
      select(1, .data$dose, .data$value) %>%
      group_by(.data$id, .data$dose) %>%
      summarize(data = list(.data$value)) %>%
      ungroup() %>%
      mutate(
        dose = as.numeric(.data$dose),
        id = as.character(.data$id)
      )
  )
}
