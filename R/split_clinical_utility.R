run_split_cu <- function(data, size) {

  decision_proportions <- dplyr::tibble(
    decision = 1:0,
    proportion_decision = c(
      mean(data$decision),
      1 - mean(data$decision)
    )
  )

  split_row_ids <- generate_split_dataset(
    data = data,
    size = size
  )

  concordant_decision <- data |>
    dplyr::filter(rowId %in% split_row_ids) |>
    dplyr::group_by(decision) |>
    tidyr::nest() |>
    dplyr::mutate(
      expected_outcome = purrr::map_dbl(
        data,
        impute_outcomes,
        decision = decision
      ),
      proportion_treated = purrr::map_dbl(
        data,
        \(x) sum(x$treatment) / nrow(x)
      ),
      treatment = decision
    ) |>
    dplyr::select(-data) |>
    dplyr::ungroup()

  discordant_decision <- data |>
    dplyr::filter(!(rowId %in% split_row_ids)) |>
    dplyr::group_by(decision) |>
    tidyr::nest() |>
    dplyr::mutate(
      expected_outcome = purrr::map_dbl(
        data,
        impute_outcomes_disc,
        decision = decision
      ),
      proportion_treated = purrr::map_dbl(
        data,
        \(x) sum(x$treatment) / nrow(x)
      ),
      treatment = as.numeric(!as.logical(decision))
    ) |>
    dplyr::select(-data) |>
    dplyr::ungroup()

  decision_utility <- concordant_decision |>
    dplyr::bind_rows(discordant_decision) |>
    dplyr::arrange(decision) |>
    dplyr::relocate(treatment, .after = decision) |>
    dplyr::relocate(proportion_treated, .after = treatment) |>
    dplyr::left_join(
      decision_proportions,
      by = "decision"
    ) |>
    dplyr::mutate(
      result = proportion_decision *
        (treatment * proportion_treated +
           (1 - treatment) * (1 - proportion_treated)) *
        expected_outcome
    ) |>
    dplyr::pull(result) |>
    sum()

  standard_of_care_utility <- mean(data$outcome)
 standard_of_care_utility - decision_utility
  
}

generate_split_dataset <- function(data, size) {

  if (size < 0 || size > 1) {
    stop("Size must be a number between 0 and 1.")
  }

  data |>
    dplyr::sample_frac(size) |>
    dplyr::pull(rowId) |>
    sort()

}

impute_outcomes <- function(data, decision) {
  data |>
    dplyr::filter(decision == treatment) |>
    dplyr::pull(outcome) |>
    mean()
}


impute_outcomes_disc <- function(data, decision) {
  estimand <- ifelse(
    decision == 0,
    yes = "ATT",
    no = "ATC"
  )

  matching <- MatchIt::matchit(
    treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = data,
    estimand = estimand,
    distance = "mahalanobis",
    ratio = 1,
    replace = TRUE
  )

  MatchIt::get_matches(matching) |>
    dplyr::filter(treatment == decision) |>
    dplyr::pull(outcome) |>
    mean()
}



# result <- furrr::future_map_dbl(1:60, ~ {
#   analysis_data |>
#     run_split_cu(size = .5)
# },
# .progress = TRUE,
# .options = furrr::furrr_options(seed = TRUE)
# )
# 
# message("\n")
# mean(result)
