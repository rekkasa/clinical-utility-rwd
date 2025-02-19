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

match_with_leftout <- function(data, analysis_data, leftout_row_ids) {

  analysis_data |>
    dplyr::filter(
      rowId %in% leftout_row_ids,
      treatment == 1
    ) |>
    dplyr::bind_rows(data) |>
    dplyr::mutate(propensityScore = plogis(fitted_ps)) |>
    CohortMethod::matchOnPs() |> dim()
  
}

# ---- DELETE ---->

data <- decision_treat_boot |>
  dplyr::left_join(analysis_data, by = "rowId")  |>
  dplyr::filter(treatment == 0)

# <---- DELETE ----

decision_treat_boot |>
  dplyr::left_join(analysis_data, by = "rowId")  |>
  dplyr::group_by(treatment) |>
  tidyr::nest() |>
  dplyr::mutate(
    data = ifelse(
      treatment == 0,
      purrr::map(
        data,
        ~ match_with_leftout(.x)
      ),
      data
    )
  )


# Find matches on PS


# Draw bootstrap sample for decision:no_treat


# Find matches on PS


# Calculate clinical utility


# Calculate true clinical utility


# Compare
