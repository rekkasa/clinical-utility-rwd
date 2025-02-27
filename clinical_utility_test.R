source("test.R")
source("helper.R")
future::plan(future::multisession, workers = 15)

cutoff_value <- .05

# Estimate propensity scores and risk
settings$databaseSettings$numberOfObservations <- 2e4

analysis_data <- SimulateHte::runDataGeneration(
  databaseSettings = settings$databaseSettings,
  propensitySettings = settings$propensitySettings,
  baselineRiskSettings = settings$baselineRiskSettings,
  treatmentEffectSettings = settings$treatmentEffectSettings
) |>
  dplyr::select(
    rowId,
    paste0("x", 1:10),
    treatment, outcome,
    outcomeTreated, outcomeUntreated
  )

ps_model <- glm(
  treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
  family = "binomial",
  data = analysis_data
)

risk_model <- glm(
  outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
  family = "binomial",
  data = analysis_data
)

analysis_data <- analysis_data |>
  dplyr::mutate(
    fitted_risk = risk_model$linear.predictors,
    decision = as.numeric(plogis(fitted_risk) > cutoff_value)
  )

# Re-estimate PS for each decision group
analysis_data <- analysis_data |>
  dplyr::group_by(decision) |>
  tidyr::nest() |>
  dplyr::mutate(
    data = purrr::map(
      data,
      ~ dplyr::mutate(
        .x,
        fitted_ps = glm(
          treatment ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
          family = "binomial",
          data = .x
        )$linear.predictors
      )
    )
  ) |>
  dplyr::ungroup(decision) |>
  tidyr::unnest(data)

clinical_utility_treat <- analysis_data |>
  dplyr::filter(decision == 1) |>
  compute_sth(n_boot = 200)

clinical_utility_no_treat <- analysis_data |>
  dplyr::filter(decision == 0) |>
  compute_sth(n_boot = 200)

threshold_clinical_utility <-
  mean(analysis_data$decision) * clinical_utility_treat +
  (1 - mean(analysis_data$decision)) * clinical_utility_no_treat

threshold_cu <- mean(analysis_data$outcome) - threshold_clinical_utility

message(
  "Threshold CU: ",
  threshold_cu
)

population <- readRDS("data/raw/population.rds")

population <- population |>
  dplyr::mutate(
    fitted_risk = predict(
      risk_model,
      newdata = population,
      type = "response"
    ),
    decision = ifelse(fitted_risk > cutoff_value, 1, 0)
  )

true_threshold_clinical_utility <- mean(
  population$decision * population$outcomeTreated +
    (1 - population$decision) * population$outcomeUntreated
)

true_threshold_cu <- mean(population$outcome) - true_threshold_clinical_utility
message(
  "True threshold CU: ",
  true_threshold_cu
)

message(
  "Deviation from the truth: ",
  paste0(
    100 * round(
    (abs(threshold_cu - true_threshold_cu) / true_threshold_cu),
    digits = 4
    ),
    "%"
  )
)
