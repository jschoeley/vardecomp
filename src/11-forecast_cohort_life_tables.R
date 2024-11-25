# Forecast cohort life tables
#
# For birth cohorts 1900 through 2019 we forecast their cohort
# life tables based on the diagonal of the associated period life
# tables. The forecasts are stochastic and based on the 30 years of data
# prior to the birth of a cohort.

# Init ------------------------------------------------------------

library(yaml)
library(readr)
library(dplyr)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
#  config = 'cfg/config.yaml',
  global = './src/00-global_objects.R',
  plc = './src/00-penalized_lee_carter_functions.R',
  lcfit_folder = './out/10-lcfit'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  lcforecastperiod = './out/11-lcforecastperiod.rds',
  lcforecastcohort = './out/11-lcforecastcohort.rds'
)

# global configuration
#config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- within(list(), {
  age_start = 0
  age_end = 100
  nsim = 100
  seed = 1987
  # sliding window fitting period
  years = 1861:1900
  increments = 0:(2022-max(years))
  prob_crisis = c(0.0362, 0.3846)
  magnitude_crisis = c(
    24, 34, 46, 57, 10, 13, 8, 7, 13, 20, 31, 6, 8
  )
})

# list containers for analysis artifacts
dat <- list()

# Functions -------------------------------------------------------

# global objects and functions
source(paths$input$global)
source(paths$input$plc)

# Load data -------------------------------------------------------

# lee carter model fits
LC_fit <- lapply(list.files(paths$input$lcfit_folder, full.names = TRUE), readRDS)
LC_fit <- unlist(LC_fit, recursive = FALSE)

# Period Forecast -------------------------------------------------

# perform a 100 year period forecast starting at each cohort birth year
# this allows to extract a 100 year cohort diagonal from the period*age
# surface

LC_forecast_period <- LC_fit

LC_forecast_period <- lapply(LC_forecast_period, function (l) {

  cat(l$meta$population, '\n')

  fcst <- PLCforecast(
    theta = l$model_parameters, h = length(l$model_parameters$ax),
    nsim = cnst$nsim,
    sd_estimation = 'classic',
    drift_estimation = 'classic',
    jumpoff_estimation = 'robust',
    kt_exclude = l$meta$kt_exclude,
    p_crisis = cnst$prob_crisis,
    m_crisis = cnst$magnitude_crisis
  )

  return(
    list(
      Eta_forecast_sim = fcst$Eta_forecast_sim,
      kt_drift = fcst$kt_drift,
      dkt_sd = fcst$dkt_sd,
      meta = l$meta
    )
  )

})

# demonstrate forecast with crises
LC_forecast_period[[1]]$dkt_sd
PlotMatrix(LC_forecast_period[[11]]$Eta_forecast_sim[,,1],
           type = 'c', N = 10)
plot(exp(LC_forecast_period$switzerland_2019$Eta_forecast_sim[90,,1]))

# Cohort Forecast -------------------------------------------------

# convert period to cohort forecasts by extracting the cohort diagonal

LC_forecast_cohort <- LC_forecast_period

LC_forecast_cohort <- lapply(LC_forecast_period, function (l) {

  cat(l$meta$population, '\n')

  # extract cohort
  Eta_forecast_sim <-
    apply(l$Eta_forecast_sim, 3, function (X) { diag(X) })

  return(
    list(
      cohort = l$meta$years |> as.numeric() |> max() + 1,
      Eta_forecast_sim = Eta_forecast_sim,
      kt_drift = l$kt_drift,
      dkt_sd = l$dkt_sd,
      meta = l$meta
    )
  )
})

# Export ----------------------------------------------------------

saveRDS(LC_forecast_period, paths$output$lcforecastperiod)
saveRDS(LC_forecast_cohort, paths$output$lcforecastcohort)
