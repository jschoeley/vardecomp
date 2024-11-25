# Implement Penalized Lee-Carter Poisson Mortality Forecast
#
# Forecast age specific mortality rates for using a 30 year
# sliding window fitting period starting 1871 and ending 2019.
# These fits will be used to derive cohort mortality forecasts
# using the data available in the 30 years prior to the birth
# of a cohort. Neighboring Lee-Carter fits are penalized to be
# similar.

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
  config = './cfg/config.yaml',
  global = './src/00-global_objects.R',
  plc = './src/00-penalized_lee_carter_functions.R',
  DxEx = './dat/hmd'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  data = './dat/output_data.rds',
  lcfit = './out/10-lcfit'
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
  years = 1871:1900
  increments = 0:(2019-max(years))
})

# list containers for analysis artifacts
dat <- list()

# Functions -------------------------------------------------------

# global objects and functions
source(paths$input$global)
source(paths$input$plc)

# Load data -------------------------------------------------------

# prepare list of country specific data frames with deaths and exposures
# by age and year
countries <- list.dirs(paths$input$DxEx, full.names = FALSE)[-1]
input_data <- lapply(countries, function (l) {
  list(
    country = l,
    path_Dxt = paste0(paths$input$DxEx, '/', l, '/Deaths_1x1.txt'),
    path_Ext = paste0(paths$input$DxEx, '/', l, '/Exposures_1x1.txt')
  )
})

input_data <- lapply(input_data, function (l) {

  Dxt <- read_table(l$path_Dxt, skip = 2, na = '.')
  Ext <- read_table(l$path_Ext, skip = 2, na = '.')

  DxEx <-
    left_join(Dxt, Ext, by = c('Year', 'Age')) |>
    select(Year, Age, Dx = Male.x, Ex = Male.y) |>
    mutate(Age = as.integer(Age),
           Age = ifelse(is.na(Age), 110L, Age),
           Dx = floor(Dx)) |>
    filter(Age %in% cnst$age_start:cnst$age_end)

  list(
    country = l$country,
    DxEx = DxEx
  )
})

names(input_data) <- countries

# Fit LC ----------------------------------------------------------

# define fitting periods. 30 year sliding window with
# single year increments. First fitting window: 1871:1900. Last
# fitting window: 1990:2019
fitting_periods <- vector('list', length(cnst$increments))
for (i in 1:length(fitting_periods)) {
  fitting_periods[[i]][['years']] <- cnst$years + cnst$increments[i]
}

LC_fit <- vector('list', length(countries)*length(fitting_periods))
country_period_names <- expand.grid(
  sapply(fitting_periods, function (l) max(l$years)), countries
) |>
  mutate(name = paste(Var2, Var1, sep = '_')) |>
  pull(name)
names(LC_fit) <- country_period_names

# configure fitting process
LC_fit_config <- PLCfitConfig(
  # ensure smoothness of ax over age
  lambda_ax = 1e-6,
  # ensure smoothness of bx over age
  lambda_bx = 1e8,
  lambda_kt = 0,
  lambda_ridge = 1e-6,
  # penalties for deviations of parameters between
  # neighboring years
  lambda_ax_target = 1e2,
  lambda_bx_target = 1e2,
  lambda_kt_target = 1e2,
  dev_stop_crit = 1e-3
)

for (country in 1:length(input_data)) {

  country_name <- input_data[[country]][['country']]

  for (period in 1:length(fitting_periods)) {

    fitting_years <- fitting_periods[[period]][['years']]
    population <- paste(country_name, max(fitting_years), sep = '_')

    cat('\n', country_name, 'Fit: ', min(fitting_years), ':', max(fitting_years), '\n')

    ages <- cnst$age_start:cnst$age_end
    df <- input_data[[country]][['DxEx']] |>
      filter(Year %in% fitting_years)

    N = length(ages)
    m = length(fitting_years)
    Dxt <- matrix(df$Dx, nrow = N, ncol = m,
                  dimnames = list(Age = ages, Year = fitting_years))
    Ext <- matrix(df$Ex, nrow = N, ncol = m,
                  dimnames = list(Age = ages, Year = fitting_years))

    year_exclusion <- which(fitting_years %in% c(1914:1921, 1939:1945))
    W <- PLCcreateWeightMatrix(Dxt, Ext, nacols = year_exclusion)

    LC_fit[[population]] <- tryCatch({
      if (period > 1) {
        # use fitted values from previous years as starting parameters
        LC_fit_config$init_pars <-
          LC_fit[[paste0(country_name, '_', max(fitting_years-1))]][[
            'unconstrained_model_parameters'
          ]]
        # shift kt 1 step forward as we are now fitting 1 year later
        LC_fit_config$init_pars$kt <- c(LC_fit_config$init_pars$kt[-1],0)
        # penalize differences from previous fit
        LC_fit_config$ax_target <- LC_fit_config$init_pars$ax
        LC_fit_config$bx_target <- LC_fit_config$init_pars$bx
        LC_fit_config$kt_target <- LC_fit_config$init_pars$kt
      }
      PLCfit(
        Dxt = Dxt, Ext = Ext, W = W,
        label = paste(country_name, 'Male', min(fitting_years),
                      max(fitting_years), sep = '-'),
        config = LC_fit_config
      )}, error = function (error) {NULL}
    )
  }

  # save fit results for single country
  saveRDS(LC_fit[grepl(country_name, names(LC_fit))],
          file = paste0(paths$output$lcfit, '/', country_name, '.rds'))

}
