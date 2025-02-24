#!/usr/bin/env Rscript

library(SimulateHte)

databaseSettings <- createDatabaseSettings(
  numberOfObservations = 1e6,
  numberOfCovariates = 10,
  covariateDistributionSettings = list(
    SimulateHte::createMultivariateNormalDistributionSettings(
      mean = rep(0, 10),
      covariance = matrix(
        c(
          c(1, 0, 0, 0, .2, rep(0, 5)),
          c(0, 1, 0, 0, 0, .9, rep(0, 4)),
          c(0, 0, 1, rep(0, 4), .2, 0, 0),
          c(0, 0, 0, 1, rep(0, 4), .9, 0),
          c(.2, 0, 0, 0, 1, rep(0, 5)),
          c(0, .9, 0, 0, 0, 1, rep(0, 4)),
          c(rep(0, 6), 1, rep(0, 3)),
          c(0, 0, .2, rep(0, 4), 1, 0, 0),
          c(rep(0, 3), .9, rep(0, 4), 1, 0),
          c(rep(0, 9), 1)
        ),
        ncol = 10
      )
    )
  )
)

propensitySettings <- createPropensitySettings(
  type = "binary",
  modelSettings = createModelSettings(
    constant = -1.897,
    modelMatrix = diag(10),
    transformationSettings = replicate(n = 10, identity, simplify = FALSE),
    coefficients = c(
      .8, -.25, .6, -.4,
      -.8, -.5, .7,
      rep(0, 3)
    )
  )
)

baselineRiskSettings <- SimulateHte::createBaselineRiskSettings(
  type = "binary",
  modelSettings = SimulateHte::createModelSettings(
    constant = -2,
    modelMatrix = diag(10),
    transformationSettings = replicate(n = 10, identity, simplify = FALSE),
    coefficients = c(
      .3, -.36, -.73, -.2,
      0, 0, 0,
      .71, -.19, .26
    )
  )
)

settings <- list(
  databaseSettings = databaseSettings,
  propensitySettings = propensitySettings,
  baselineRiskSettings = baselineRiskSettings
)

saveRDS(settings, "data/processed/leacySettings.rds")
