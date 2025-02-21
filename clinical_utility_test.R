source("test.R")
source("helper.R")
future::plan(future::multisession, workers = 12)

cutoff_value <- .1

# Estimate propensity scores and risk

analysis_data <- simulated_data |>
  dplyr::sample_n(1e5) |>
  dplyr::select(
    rowId,
    paste0("x", 1:8),
    treatment, outcome,
    outcomeTreated, outcomeUntreated
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

pp <- analysis_data |>
  dplyr::filter(decision == 0) |>
  compute_sth(n_boot = 500)

ll <- simulated_data |>
  dplyr::mutate(
    fitted_risk = predict(
      risk_model,
      newdata = simulated_data,
      type = "response"
    ),
    decision = as.numeric(fitted_risk > cutoff_value)
  ) |>
  dplyr::filter(decision == 0) |>
  dplyr::pull(outcomeUntreated) |>
  mean()


pp - ll
