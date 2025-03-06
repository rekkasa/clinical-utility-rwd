source("test.R")
source("helper.R")
future::plan(future::multisession, workers = 4)

cutoff_value <- .05

# Estimate propensity scores and risk
settings$databaseSettings$numberOfObservations <- 1e4

population <- readRDS("data/raw/population.rds")
message("Read population")

n_replications <- 5
n_boot <- 10
threshold_cu_all <- threshold_cu_decision <-
  true_threshold_cu <- c(0, n_replications)

for (i in 1:n_replications) {

  message(paste("Starting analysis:", i))
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

  message("Bootstrap: all")
  message("Computing clinical utility for decision: treat")

  clinical_utility_treat <- analysis_data |>
    dplyr::filter(decision == 1) |>
    bootstrap_clinical_utility(n_boot = n_boot)


  message("\nComputing clinical utility for decision: no treat")

  clinical_utility_no_treat <- analysis_data |>
    dplyr::filter(decision == 0) |>
    bootstrap_clinical_utility(n_boot = n_boot)


  threshold_clinical_utility <-
    mean(analysis_data$decision) * clinical_utility_treat +
    (1 - mean(analysis_data$decision)) * clinical_utility_no_treat

  threshold_cu_all[i] <- mean(analysis_data$outcome) - threshold_clinical_utility

  message("\nComputed clinical utility for proposed rule")

  message("Bootstrap: decision")
  message("Computing clinical utility for decision: treat")

  clinical_utility_treat <- analysis_data |>
    dplyr::filter(decision == 1) |>
    bootstrap_clinical_utility(n_boot = n_boot, bootstrap = "decision")


  message("\nComputing clinical utility for decision: no treat")

  clinical_utility_no_treat <- analysis_data |>
    dplyr::filter(decision == 0) |>
    bootstrap_clinical_utility(n_boot = n_boot, bootstrap = "decision")


  threshold_clinical_utility <-
    mean(analysis_data$decision) * clinical_utility_treat +
    (1 - mean(analysis_data$decision)) * clinical_utility_no_treat

  threshold_cu_decision[i] <- mean(analysis_data$outcome) - threshold_clinical_utility

  message("\nComputed clinical utility for proposed rule")

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

  true_threshold_cu[i] <- mean(population$outcome) - true_threshold_clinical_utility
  message("Computed true clinical utility for proposed rule")

}


message("Saving results")
readr::write_csv(
  x = data.frame(
    boot_all = threshold_cu_all,
    boot_decision = threshold_cu_decision,
    actual = true_threshold_cu
  ),
  file = file.path(
    "data/processed",
    paste0("boot_all", ".csv")
  )
)
