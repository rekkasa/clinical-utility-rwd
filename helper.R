bootstrap_clinical_utility <- function(
  data,
  n_boot = 50,
  bootstrap = "all"
) {

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
      calculate_boot_result(bootstrap = bootstrap)
  }, .progress = TRUE)

  return(mean(res))

}

calculate_boot_result <- function(
  data,
  bootstrap = "all",
  ...
) {

  if (bootstrap == "all") {
    boot_sample_size <- nrow(data)

    decision_boot <- data.frame(
      rowId = sample(
        x = data |> dplyr::pull(rowId),
        size = boot_sample_size,
        replace = TRUE
      )
    )

    sample_to_match <- decision_boot
  } else if (bootstrap == "decision") {
    boot_sample_size <- sum(data$treatment != data$decision)

    decision_boot <- data.frame(
      rowId = sample(
        x = data |>
          dplyr::filter(decision != treatment) |>
          dplyr::pull(rowId),
        size = boot_sample_size,
        replace = TRUE
      )
    )

    sample_to_match <- data |>
      dplyr::filter(treatment == decision) |>
      dplyr::select(rowId) |>
      dplyr::bind_rows(decision_boot)
  } else if (bootstrap == "leftout") {

    sample_to_match <- decision_boot
  }

  estimand <- ifelse(unique(data$decision) == 0, "ATT", "ATC")

  matching <- MatchIt::matchit(
    treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = sample_to_match |> dplyr::left_join(data, by = "rowId"),
    estimand = estimand,
    distance = "mahalanobis",
    ratio = 1,
    replace = TRUE,
    # caliper = .1,
    ...
  )

  matched_data <- MatchIt::get_matches(matching)

  mean_outcome_disc <- matched_data |>
    dplyr::filter(treatment == decision) |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(mean_outcome = mean(outcome)) |>
    dplyr::mutate(treatment = as.numeric(!treatment))

  mean_outcome_conc <- sample_to_match |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::filter(treatment == decision) |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(mean_outcome = mean(outcome))

  result <- sample_to_match |>
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

calculate_boot_result_leftout <- function(
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

  disc_row_ids <- decision_boot |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::filter(treatment == disc_treatment) |>
    dplyr::pull(rowId)

  leftout_row_ids <- setdiff(
    data$rowId,
    unique(decision_boot$rowId)
  )

  matched_row_ids <- decision_boot |>
    dplyr::left_join(
      data,
      by = "rowId"
    ) |>
    dplyr::filter(treatment == disc_treatment) |>
    dplyr::pull(rowId) |>
    match_with_leftout(
      analysis_data = data,
      leftout_row_ids = leftout_row_ids
    )

  overall_outcome_conc <- decision_boot |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::filter(treatment == decision) |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(
      average_outcome = mean(outcome),
      n_total = dplyr::n()
    )

  mean_outcome_disc <- data.frame(rowId = matched_row_ids) |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::pull(outcome) |>
    mean()

  overall_outcome <- data.frame(
    treatment = disc_treatment,
    average_outcome = mean_outcome_disc,
    n_total = length(disc_row_ids)
  ) |>
    dplyr::bind_rows(overall_outcome_conc)

  sum(
    overall_outcome$average_outcome *
      overall_outcome$n_total
  ) /
    sum(overall_outcome$n_total)
}

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

  matched_data <- MatchIt::get_matches(matching)

  matched_data |>
    dplyr::filter(treatment == decision) |>
    dplyr::pull(rowId)
}
