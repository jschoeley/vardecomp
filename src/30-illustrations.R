# Survivorship illustrations
#
# How does the script do it?

# Init ------------------------------------------------------------

library(yaml)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Constants -------------------------------------------------------

# input and output paths
setwd('.')
paths <- list()
paths$input <- list(
  tmpdir = './tmp',
  config = './cfg/config.yaml',
  global = './src/00-global_objects.R',
  lcforecastcohort = './out/11-lcforecastcohort.rds'
)
paths$output <- list(
  tmpdir = paths$input$tmpdir,
  fig = './out/'
)

# global configuration
# config <- read_yaml(paths$input$config)

# constants specific to this analysis
cnst <- within(list(), {
  seed = 1987
  nsim = 1e4
})

set.seed(cnst$seed)

# list containers for analysis artifacts
dat <- list()

# Load data -------------------------------------------------------

cohort_forecasts <- readRDS(paths$input$lcforecastcohort)

# Functions -------------------------------------------------------

# global objects and functions
source(paths$input$global)

# sample from siler distribution
rsiler <- function (
    n = 100, a1 = 0.001, b1 = -0.2, a2 = 0.0003, a3 = 0.0001, b3 = 0.14
) {
  require(flexsurv)
  X1 <- rgompertz(n = n, rate = a1, shape = b1)
  X2 <- rexp(n = n, rate = a2)
  X3 <- rgompertz(n = n, rate = a3, shape = b3)
  X <- pmin(X1, X2, X3)
  return(X)
}

# density, survival, and hazard of siler distribution
fsiler <- function (
    x, a1 = 0.001, b1 = -0.2, a2 = 0.0003, a3 = 0.0001, b3 = 0.14
) {
  require(flexsurv)
  # survival function of Siler is S1(x)*S2(x)*S3(x)
  S1 <- pgompertz(x, rate = a1, shape = b1, lower.tail = FALSE)
  S2 <- pexp(x, rate = a2, lower.tail = FALSE)
  S3 <- pgompertz(x, rate = a3, shape = b3, lower.tail = FALSE)
  S <- S1*S2*S3
  # Siler hazard
  h <- a1*exp(b1*x) + a2 + a3*exp(b3*x)
  # Siler density
  d <- h*S
  
  return(list(density = d, hazard = h, survival = S))
}

# Siler simulations -----------------------------------------------

sim <- list(
  outsurvival = NULL,
  extremes = NULL,
  forecasts = NULL
)

sim$outsurvival <- tibble(
  a = rsiler(n = 1e5, a3 = 1e-6),
  b = rsiler(n = 1e5, a3 = 4e-6)
) |>
  pivot_longer(everything(), names_to = 'id')
sim$extremes <- tibble(
  no_aging = rsiler(n = 1e4, b3 = 0),
  famine = rsiler(n = 1e4, a1 = 0.1)
) |>
  pivot_longer(everything(), names_to = 'id')
sim$forecasts <- tibble(
  a3 = exp(rnorm(1e2, log(1e-6), 1)), b3 = rnorm(1e2, 0.14, sd = 0.01)
) |>
  mutate(id = as.character(1:1e2)) |>
  group_by(id) |>
  group_modify(~{
    tibble(
      value = rsiler(n = 1e3, a3 = .x$a3, b3 = .x$b3)
    )
  }) |>
  ungroup()

sim$all <- bind_rows(
  forecasts = sim$forecasts,
  extremes = sim$extremes,
  outsurvival = sim$outsurvival,
  .id = 'stratum'
) |>
  arrange(stratum, id, value)

# Lee Carter simualtions ------------------------------------------

ew2019 <- cohort_forecasts$englandwales_2019$Eta_forecast_sim
dimnames(ew2019) <- list(age = 1:101, sim = 1:100)

ew2019 <-
  ew2019 |>
  apply(MARGIN = 2, function (logmx) {
    # survival lx
    exp(cumsum(-exp(logmx)))
  }) |>
  as.data.frame.table() |>
  as_tibble() |>
  mutate(
    stratum = 'ew2019',
    sim = as.character(sim),
    age = as.double(as.character(age))) |> 
  select(stratum, id = sim, value = age, p = Freq)


# Survival curves -------------------------------------------------

survival <-
  sim$all |>
  group_by(stratum, id) |>
  mutate(
    # observed survival (KM)
    p = 1 - 1:n()/n()
  ) |>
  ungroup() |>
  bind_rows(ew2019)

survivorship <-
  survival |>
  ggplot(aes(x = value, y = p)) +
  geom_step(
    aes(group = id), color = 'black', alpha = 0.5,
    data = . %>% filter(stratum == 'forecasts')
  ) +
  geom_line(
    aes(group = id), color = 'blue', alpha = 0.5,
    data = . %>% filter(stratum == 'ew2019')
  ) +
  geom_step(
    aes(group = id), color = 'red', alpha = 0.5,
    data = . %>% filter(stratum == 'extremes')
  ) +
  geom_step(
    aes(group = id), color = 'yellow',
    data = . %>% filter(stratum == 'outsurvival')
  ) +
  coord_cartesian(xlim = c(0, 100), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(x = 'Age', y = 'Survivorship') +
  MyGGplotTheme()


# Statistics of interest ------------------------------------------

summaries <-
  sim$all |>
  mutate(a_survives_b = a > b) |>
  summarise(
    e0_a = mean(a),
    e0_b = mean(b),
    p_a_survives_b = sum(a_survives_b) / n()
  )

# within group variability higher than between group
summaries <-
  sim$outsurvival |>
  pivot_wider(names_from = id, values_from = value) |>
  unnest(c(a, b)) |>
  mutate(a_survives_b = a > b) |>
  summarise(
    e0_a = mean(a),
    e0_b = mean(b),
    p_a_survives_b = sum(a_survives_b) / n(),
    mad_A = mean(abs(
      sample(a, n(), replace = TRUE) -
        sample(a, n(), replace = TRUE)
    )),
    mad_B = mean(abs(
      sample(b, n(), replace = TRUE) -
        sample(b, n(), replace = TRUE)
    ))
  )

# Export ----------------------------------------------------------

ExportFigure(
  survivorship, path = paths$output$fig,
  filename = '30-survivorship.rds',
  device = 'pdf',
  width = 170, height = 130
)
