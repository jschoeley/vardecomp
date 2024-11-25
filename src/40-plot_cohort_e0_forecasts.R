library(ggplot2)
library(dplyr)

# Cohort e0: Observed vs. prediction error ------------------------

LC_forecast_cohort <- readRDS('./out/11-lcforecastcohort.rds')

source('src/00-global_objects.R')
source('src/00-penalized_lee_carter_functions.R')

predicted_cohort_ex <-
  lapply(LC_forecast_cohort, function (l) {
    ex_pred <- apply(exp(l$Eta_forecast_sim), 2,
                     LifeExpectancyFromMortality)
    colnames(ex_pred) <- 1:ncol(l$Eta_forecast_sim)
    rownames(ex_pred) <- l$meta$ages
    ex_pred <- as.data.frame(as.table(ex_pred))
    colnames(ex_pred) <- c('age', 'sim', 'ex')
    ex_pred$age <- as.integer(ex_pred$age)
    ex_pred$sim <- as.integer(ex_pred$sim)

    data.frame(cohort = l$cohort, ex_pred)
  })

bind_rows(predicted_cohort_ex, .id = 'country') |>
  mutate(country = gsub(pattern = '_.+', replacement = '', x = country)) |>
  filter(
    age == 1
  ) |>
  group_by(cohort, country) |>
  summarise(
    ex_avg = mean(ex),
    ex_lo = quantile(ex, 0.025, na.rm = T),
    ex_hi = quantile(ex, 0.975, na.rm = T)
  ) |>
  ggplot(aes(x = cohort)) +
  geom_ribbon(
    aes(ymin = ex_lo, ymax = ex_hi),
    color = NA, fill = 'grey'
  ) +
  geom_line(aes(y = ex_avg)) +
  facet_wrap(~country, nrow = 2) +
  scale_y_continuous(breaks = seq(0, 100, 5)) +
  MyGGplotTheme()
