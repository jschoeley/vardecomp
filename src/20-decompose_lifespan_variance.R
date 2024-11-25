# Lifespan Variance Decomposition
#
# How does the script do it?

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
  lcforecastperiod = './out/11-lcforecastperiod.rds',
  lcforecastcohort = './out/11-lcforecastcohort.rds',
  plc = './src/00-penalized_lee_carter_functions.R'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  vardecomp = './out/20-vardecomp.rds'
)

# global configuration
# config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- within(list(), {})

# list containers for analysis artifacts
dat <- list()

# Functions -------------------------------------------------------

# global objects and functions
source(paths$input$global)
source(paths$input$plc)

# Load data -------------------------------------------------------

LC_forecast_cohort <- readRDS(paths$input$lcforecastcohort)

# Variance Decomposition ------------------------------------------

cohort_lifespan_variance_decomp <- LC_forecast_cohort

cohort_lifespan_variance_decomp <- lapply(LC_forecast_cohort, function (l) {

  cat(l$meta$population, '\n')

  Mcx_sim <- exp(l$Eta_forecast_sim)

  # e0_i = E(X|Mx_i):
  # Cohort life expectancy conditional simulated LT
  ex_sim <- apply(Mcx_sim, 2, LifeExpectancyFromMortality)
  e0_sim <- ex_sim[1,] # at age 0

  # e0 = E(X):
  # Cohort life expectancy averaged over all simulated LTs
  ex_avg <- rowMeans(ex_sim)
  e0_avg <- ex_avg[1] # at age 0

  # var0_i = Var(X|Mx_i):
  # Cohort lifespan variance conditional on simulated LT
  varx_sim <- apply(Mcx_sim, 2, LifespanVarianceFromMortality)
  var0_sim <- varx_sim[1,] # at age 0

  # Evar0 = E_i[Var(X|Mx_i)]:
  # Average cohort lifespan variance over all simulated LTs,
  # a.k.a. within-lifetable-variance
  # a.k.a. expected cohort life-span variance
  varx_avg <- rowMeans(varx_sim)
  var0_avg <- varx_avg[1] # at age 0

  # vare0 = Var_i[E(X|Mx_i):
  # Variance of the conditional cohort life expectancies
  # a.k.a. between-lifetable-variance
  # a.k.a. forecasting variance of cohort life tables
  varex <- apply(ex_sim, 1, var)
  vare0 <- varex[1]
  
  # var_total = Var(X):
  # Total cohort lifespan variance
  var_total <- var0_avg + vare0

  # proportion of within lifetable variance on total lifespan variance
  p_var_within <- var0_avg / (var_total)*100
  cat(p_var_within, '\n')

  return(
    list(
      cohort = l$cohort,
      var_within = var0_avg,
      var_between = vare0,
      var_total = var_total,
      p_var_within = p_var_within,
      p_var_between = 100-p_var_within,
      meta = l$meta
    )
  )

})

# convert to data frame
cohort_lifespan_variance_decomp_df <- data.frame(
  country = sapply(cohort_lifespan_variance_decomp, function(l) {
    sub('-[MF].+$', replacement = '', l$meta$population)
  }),
  cohort = sapply(cohort_lifespan_variance_decomp, function(l) l$cohort),
  var_within = sapply(cohort_lifespan_variance_decomp, function(l) l$var_within),
  var_between = sapply(cohort_lifespan_variance_decomp, function(l) l$var_between),
  var_total = sapply(cohort_lifespan_variance_decomp, function(l) l$var_total),
  p_var_within = sapply(cohort_lifespan_variance_decomp, function(l) l$p_var_within),
  p_var_between = sapply(cohort_lifespan_variance_decomp, function(l) l$p_var_between)
)

# Plot Variance Decomposition -------------------------------------

fig <- list()
fig$vardecompcohort1 <- list()
fig$vardecompcohort1$data <- cohort_lifespan_variance_decomp_df

fig$vardecompcohort1$plot <-
  ggplot(fig$vardecompcohort1$data) +
  aes(x = cohort, y = p_var_between) +
  geom_line() +
  scale_x_continuous(breaks = seq(1900, 2020, 20)) +
  scale_y_continuous(breaks = seq(0, 100, 0.5), limits = c(0, NA)) +
  labs(
    title = 'The influence of forecasting uncertainty on cohort lifespan uncertainty',
    subtitle = 'Share of total cohort lifespan variance explained by variance in forecasted life-tables',
    x = 'Cohort',
    y = '%'
  ) +
  facet_wrap(~country, nrow = 2) +
  coord_cartesian(expand = FALSE) +
  MyGGplotTheme(grid = 'y', panel_border = TRUE)

fig$vardecompcohort1$plot

ExportFigure(
  fig$vardecompcohort1$plot,
  path = paths$output$tmpdir,
  '20-vardecompcohort1',
  device = 'pdf'
)

fig$vardecompcohort2$data <-
  fig$vardecompcohort1$data |>
  filter(!country %in% c('denmark', 'belgium')) |>
  select(-var_total) |>
  tidyr::pivot_longer(c(var_within, var_between),
                      names_to = 'Variance component') |>
  mutate(
    `Variance component` = factor(`Variance component`,
                                  rev(c('var_between', 'var_within')),
                                  rev(c('forecasting', 'lifetable')))
  )

fig$vardecompcohort2$plot <-
  ggplot(fig$vardecompcohort2$data) +
  aes(x = cohort) +
  geom_area(aes(y = value, fill = `Variance component`)) +
  scale_x_continuous(breaks = seq(1900, 2020, 20)) +
  scale_y_continuous(
    breaks = seq(0, 100, 5)^2, limits = c(0, NA),
    sec.axis = sec_axis('sqrt', breaks = seq(0,100,5), name = 'Standard deviation')
  ) +
  scale_fill_manual(values = c('grey70', '#D7000C')) +
  labs(
    title = 'The influence of forecasting uncertainty on cohort lifespan uncertainty',
    subtitle = 'Variance decomposition of total cohort lifespan uncertainty',
    x = 'Cohort',
    y = 'Variance'
  ) +
  facet_wrap(~country, nrow = 2) +
  coord_cartesian(expand = FALSE) +
  MyGGplotTheme(grid = 'xy', panel_border = TRUE) +
  theme(legend.position = c(0.9, 0.9))

fig$vardecompcohort2$plot

ExportFigure(
  fig$vardecompcohort2$plot,
  path = paths$output$tmpdir,
  '20-vardecompcohort2',
  device = 'pdf'
)

saveRDS(cohort_lifespan_variance_decomp, paths$output$vardecomp)