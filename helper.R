match_with_leftout <- function(
    disc_row_ids,
    analysis_data,
    leftout_row_ids,
    ...) {
  disc_decision <- analysis_data |>
    dplyr::filter(rowId %in% disc_row_ids) |>
    dplyr::pull(treatment) |>
    unique()

  if (length(disc_decision) > 1) stop("More than one decisions in data.")
  conc_decision <- (!disc_decision) * 1

  estimand <- ifelse(disc_decision == 1, "ATT", "ATC")

  # combined_data <- analysis_data |>
  #   dplyr::filter(
  #     rowId %in% leftout_row_ids,
  #     treatment == decision
  #   ) |>
  #   dplyr::bind_rows(
  #     analysis_data |>
  #       dplyr::filter(rowId %in% disc_row_ids)
  #   )

  combined_data <- analysis_data |>
    dplyr::filter(
      rowId %in% leftout_row_ids,
      treatment == decision
    ) |>
    dplyr::bind_rows(
      data.frame(
        rowId = disc_row_ids
      ) |>
        dplyr::left_join(analysis_data, by = "rowId")
    )

  matching <- MatchIt::matchit(
    # treatment ~ fitted_ps,
    treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = combined_data,
    estimand = estimand,
    distance = "glm",
    ratio = 1,
    replace = TRUE,
    caliper = .1,
    ...
  )

  # matched_data <- MatchIt::match.data(matching)
  matched_data <- MatchIt::get_matches(matching)


  # matched_data |>
  #   dplyr::mutate(
  #     tt = ifelse(treatment == 1, "treated", "untreated")
  #   ) |>
  #   tidyr::pivot_wider(
  #     names_from = tt,
  #     values_from = rowId,
  #     id_cols = subclass
  #   ) |> 
  #   dplyr::select(treated_rowId, untreated_rowId)

  matched_data |>
    dplyr::filter(treatment == decision) |>
    dplyr::pull(rowId)
}

# data: Should only contain columns `rowId, `treatment`,
# `fitted_ps`, `decision`
compute_sth <- function(data, n_boot = 50) {
  decision <- unique(data$decision)
  if (length(decision) > 1) {
    stop("More than one decisions")
  }
  conc_treatment <- decision
  disc_treatment <- (!decision) * 1

  n_decision <- nrow(data)

  res <- rep(0, n_boot)

  res <- furrr::future_map_dbl(1:n_boot, ~ {
    data |>
      calculate_boot_result(
        n_decision = n_decision,
        disc_treatment = disc_treatment,
        conc_treatment = conc_treatment
      )
  }, .progress = TRUE)

  return(mean(res))

}

calculate_boot_result <- function(
  data,
  n_decision,
  disc_treatment,
  conc_treatment
) {

  decision_boot <- data.frame(
    rowId = sample(
      x = data |> dplyr::pull(rowId),
      size = n_decision,
      replace = TRUE
    )
  )

  estimand <- ifelse(disc_treatment == 1, "ATT", "ATC")
  matching <- MatchIt::matchit(
    treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = decision_boot |> dplyr::left_join(data, by = "rowId"),
    estimand = estimand,
    distance = "mahalanobis",
    ratio = 1,
    replace = TRUE
    # caliper = .1
  )

  matched_data <- MatchIt::get_matches(matching)

  mean_outcome_disc <- matched_data |>
    dplyr::filter(treatment == decision) |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(mean_outcome = mean(outcome)) |>
    dplyr::mutate(treatment = as.numeric(!treatment))

  mean_outcome_conc <- decision_boot |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::filter(treatment == decision) |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(mean_outcome = mean(outcome))

  result <- decision_boot |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(n_total = dplyr::n()) |>
    dplyr::left_join(
      mean_outcome_conc |>
        dplyr::bind_rows(mean_outcome_disc),
      by = "treatment"
    )

  sum(result$n_total * result$mean_outcome) / sum(result$n_total)

}
