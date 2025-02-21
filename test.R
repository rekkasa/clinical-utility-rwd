#!/usr/bin/env Rscript
library(SimulateHte)
analysisIds <- readr::read_csv(
  "https://raw.githubusercontent.com/mi-erasmusmc/HteSimulationRCT/main/data/processed/analysisIds.csv",
  col_types =  readr::cols(
    .default = readr::col_double(),
    base = readr::col_character(),
    type = readr::col_character(),
    harm = readr::col_character()
  )
)

selectedScenario <- 217    # This can be any of the possible 648 scenarios
idSettings <- analysisIds |> 
  dplyr::filter(scenario == selectedScenario)

databaseSettings <- createDatabaseSettings(
  numberOfObservations = 1e6,
  numberOfCovariates = 8,
  covariateDistributionSettings = list(
    createNormalDistributionSettings(),
    createNormalDistributionSettings(),
    createNormalDistributionSettings(),
    createNormalDistributionSettings(),
    createBinomialDistributionSettings(prob = .2),
    createBinomialDistributionSettings(prob = .2),
    createBinomialDistributionSettings(prob = .2),
    createBinomialDistributionSettings(prob = .2)
  )
)

baselineRiskSettings <- createBaselineRiskSettings(
  type = "binary",
  modelSettings = createModelSettings(
    constant = idSettings$b0,
    modelMatrix = diag(8),
    transformationSettings = replicate(n = 8, identity, simplify = FALSE),
    coefficients = idSettings |> dplyr::select(paste0("b", 1:8)) |> unlist()
  )
)

propensitySettings <- createPropensitySettings(
  type = "binary",
  modelSettings = createModelSettings(
    constant = 0,
    modelMatrix = diag(0),
    transformationSettings = list(),
    coefficients = c()
  )
)


# Functions for generating linear, quadratic, and non-monotonic deviations
# from constant relative treatment effect
createF1 <- function(c) function(x) x - c
createF2 <- function(c) function(x) (x - c)^2

# createModelSettings defines the true treatment effect function. The settings
# below define the function:
#        lp1 = g0 + g1 * (lp0 - c) + g2 * (lp0 - c)^2
# where lp1 is the true linear predictor in the treatment arm and
# lp0 is the true linear predictor in the control arm (see paper)

treatmentEffectSettings <- createTreatmentEffectSettings(
  type = "lp",
  harm = 0,
  modelSettings = createModelSettings(
    constant = idSettings$g0,
    modelMatrix = matrix(c(1, 1)),
    transformationSettings = list(
      createF1(idSettings$c),
      createF2(idSettings$c)
    ),
    coefficients = c(
      idSettings$g1,
      idSettings$g2
    )
  )
)

simulated_data <- runDataGeneration(
  databaseSettings = databaseSettings,
  baselineRiskSettings = baselineRiskSettings,
  propensitySettings = propensitySettings,
  treatmentEffectSettings = treatmentEffectSettings
)

message(
  "Difference in outcomes: ",
  mean(simulated_data$outcomeUntreated) - mean(simulated_data$outcomeTreated)
)
message(
  "Average true benefit: ",
  mean(simulated_data$trueBenefit)
)
