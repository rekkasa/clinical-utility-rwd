source("test.R")

cutoff_value <- .1

# Estimate propensity scores and risk

analysis_data <- simulated_data |>
  dplyr::select(
    rowId,
    paste0("x", 1:8),
    treatment, outcome
  )

ps_model <- glm(
  treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
  family = "binomial",
  data = analysis_data
)

risk_model <- glm(
  outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8,
  family = "binomial",
  data = analysis_data
)

analysis_data <- analysis_data |>
  dplyr::mutate(
    fitted_risk = risk_model$linear.predictors,
    fitted_ps = ps_model$linear.predictors,
    decision = as.numeric(plogis(fitted_risk) > cutoff_value)
  )

# Split in decision:treat / decision:no_treat

decision_treat_row_ids <- analysis_data |>
  dplyr::filter(decision == 1) |>
  dplyr::pull(rowId)

decision_no_treat_row_ids <- analysis_data |>
  dplyr::filter(decision == 0) |>
  dplyr::pull(rowId)


# Draw bootstrap sample for decision:treat

n_decision_treat <- length(decision_treat_row_ids)
decision_treat_boot <- data.frame(
  rowId = sample(
    x = decision_treat_row_ids,
    size = n_decision_treat,
    replace = TRUE
  )
)

match_with_leftout <- function(
  disc_row_ids,
  analysis_data,
  leftout_row_ids,
  ...
) {

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
      analysis_data |>
        dplyr::filter(rowId %in% disc_row_ids)
    )

  matching <- MatchIt::matchit(
    treatment ~ fitted_ps,
    data = combined_data,
    estimand = estimand,
    ...
  )
  matched_data <- MatchIt::match.data(matching)

  matched_data |>
    dplyr::filter(treatment == conc_decision) |>
    dplyr::pull(rowId)
}

# ---- DELETE ---->

disc_row_ids <- decision_treat_boot |>
  dplyr::left_join(analysis_data, by = "rowId") |>
  dplyr::filter(treatment == 0) |>
  dplyr::pull(rowId)

leftout_row_ids <- setdiff(
  decision_treat_row_ids,
  unique(decision_treat_boot$rowId)
)

# <---- DELETE ----

matched_row_ids <- decision_treat_boot |>
  dplyr::left_join(
    analysis_data,
    by = "rowId"
  ) |>
  dplyr::filter(treatment == 0) |>
  dplyr::pull(rowId) |>
  match_with_leftout(
    analysis_data = analysis_data,
    leftout_row_ids = leftout_row_ids
  )

overall_outcome_conc <- decision_treat_boot |>
  dplyr::filter(!(rowId %in% disc_row_ids)) |>
  dplyr::left_join(analysis_data, by = "rowId") |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(average_outcome = mean(outcome), n_total = dplyr::n())

mean_outcome_disc <- data.frame(rowId = matched_row_ids) |>
  dplyr::left_join(analysis_data) |>
  pull(outcome) |>
  mean()

overall_outcome <- data.frame(
  treatment = 0,
  average_outcome = mean_outcome_disc,
  n_total = length(disc_row_ids)
) |>
  dplyr::bind_rows(overall_outcome_conc)

sum(overall_outcome$average_outcome * overall_outcome$n_total / sum(overall_outcome$n_total))
