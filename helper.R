bootstrap_clinical_utility <- function(
  data,
  n_boot = 50,
  bootstrap = "all",
  rule_decision
) {

  res <- rep(0, n_boot)

  res <- furrr::future_map_dbl(1:n_boot, ~ {
    data |>
      calculate_boot_result(
        bootstrap = bootstrap,
        rule_decision = rule_decision
      )
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
  )

  return(mean(res))

}

calculate_boot_result <- function(
  data,
  bootstrap = "all",
  rule_decision,
  ...
) {

  if (bootstrap == "all") {

    data <- data |>
      dplyr::filter(decision == rule_decision)

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

    data <- data |>
      dplyr::filter(decision == rule_decision)

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

    result <- calculate_boot_result_leftout(
      data = data,
      rule_decision = rule_decision
    )

    return(result)

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
  rule_decision
) {

  n_decision = sum(data$decision == rule_decision)

  decision_boot <- data.frame(
    rowId = sample(
      x = data |>
        dplyr::filter(decision == rule_decision) |>
        dplyr::pull(rowId),
      size = n_decision,
      replace = TRUE
    )
  )

  disc_row_ids <- decision_boot |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::filter(treatment != rule_decision) |>
    dplyr::pull(rowId)

  leftout_row_ids <- setdiff(
    data |>
      dplyr::filter(decision == rule_decision) |>
      dplyr::pull(rowId),
    unique(decision_boot$rowId)
  )

  matched_row_ids <- decision_boot |>
    dplyr::left_join(
      data,
      by = "rowId"
    ) |>
    dplyr::filter(treatment != rule_decision) |>
    dplyr::pull(rowId) |>
    match_with_leftout(
      analysis_data = data,
      leftout_row_ids = leftout_row_ids
    )

  overall_outcome_conc <- decision_boot |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::filter(treatment == rule_decision) |>
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
    treatment = as.numeric(!(as.logical(rule_decision))),
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
      treatment == conc_decision
    ) |>
    dplyr::bind_rows(
      data.frame(
        rowId = disc_row_ids
      ) |>
        dplyr::left_join(analysis_data, by = "rowId")
    )

  matching <- MatchIt::matchit(
    treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    data = combined_data,
    estimand = estimand,
    distance = "mahalanobis",
    ratio = 1,
    replace = TRUE,
    # caliper = .1,
    ...
  )

  matched_data <- MatchIt::get_matches(matching)

  matched_data |>
    dplyr::filter(treatment == decision) |>
    dplyr::pull(rowId)
}


compute_clinical_utility <- function(data) {

  ps_model <- MatchIt::matchit(
    data = data,
    formula = treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
    estimand = "ATT",
    distance = "glm",
    method = "subclass",
    subclass = 5
  )

  matched_data <- MatchIt::match_data(ps_model)

  average_treatment_effect <- matched_data |>
    dplyr::group_by(subclass) |>
    tidyr::nest() |>
    dplyr::mutate(
      outcome_model = purrr::map(
        data,
        glm,
        formula = outcome ~ treatment,
        family = "binomial"
      )
    ) |>
    dplyr::mutate(
      treatment_effect = purrr::map_dbl(
        outcome_model,
        function(x) {
          cf <- coefficients(x)
          return(cf["treatment"])
        }
      )
    ) |>
    dplyr::pull(treatment_effect) |>
    exp() |>
    mean()


  .impute_counterfactual <- function(
    data,
    actual_treatment,
    average_treatment_effect
  ) {

    if (actual_treatment == 0) {

      data |>
        dplyr::mutate(
          counterfactual = rbinom(
            n = nrow(data),
            size = 1,
            prob = plogis(fitted_risk + log(average_treatment_effect))
          )
        ) |>
        dplyr::select(rowId, counterfactual)

    } else if (actual_treatment == 1) {

      data |>
        dplyr::mutate(
          counterfactual = rbinom(
            n = nrow(data),
            size = 1,
            prob = plogis(fitted_risk)
          )
        ) |>
        dplyr::select(rowId, counterfactual)

    }
  }

  result <- data |>
    dplyr::group_by(treatment) |>
    tidyr::nest() |>
    dplyr::mutate(
      counterfactual = purrr::map(
        data,
        actual_treatment = treatment,
        .impute_counterfactual,
        average_treatment_effect = average_treatment_effect
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(counterfactual) |>
    dplyr::ungroup(treatment) |>
    dplyr::select(-treatment) |>
    dplyr::left_join(data, by = "rowId") |>
    dplyr::mutate(
      outcome_decision = as.numeric(treatment == decision) * outcome +
        as.numeric(treatment != decision) * counterfactual
    ) |>
    dplyr::group_by(decision) |>
    tidyr::nest() |>
    dplyr::mutate(
      n = purrr::map_dbl(data, nrow)
    ) |>
    dplyr::mutate(
      average_outcome = purrr::map_dbl(
        data,
        function(x) {
          x |>
            dplyr::pull(outcome_decision) |>
            mean()
        }
      )
    ) |> 
    dplyr::select(-data) |>
    dplyr::ungroup(decision)

  overall_decision = sum(result$n * result$average_outcome) / sum(result$n)
  mean(data$outcome) - overall_decision

}
