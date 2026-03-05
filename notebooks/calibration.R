
library("koalas")
library("tidyverse")
theme_set(theme_light())

KoalasV2$new()$parameters

## Calibration of acute_duration and lifespan_acute so that 25% of animals die
## and the rest (99%) progress to chronic within 1 year
model <- KoalasV2$new()
model$set_state(S = 0, Af = 100)

model$set_parameters(
  vacc_immune_duration = Inf,
  vacc_redshed_duration = Inf,
  natural_immune_duration = Inf,
  subcinical_duration = Inf,
  subclinical_recover_proportion = 0.0,
  diseased_recover_proportion = 0.0,
  birthrate = 0.0,
  acute_duration = Inf,
  lifespan_natural = Inf,
  lifespan_acute = Inf,
  lifespan_chronic = Inf,
  relative_fecundity = 0.0,
  sensitivity = 1.0,
  specificity = 1.0,
  cure_prob_N = 1.0,
  cure_prob_I = 1.0,
  cure_prob_A = 1.0,
  cure_prob_C = 1.0,
  vaccine_efficacy = 1.0,
  vaccine_booster = 1.0,
  passive_intervention_rate = 0.0
)

model$set_parameters(
  acute_duration = 0.5,
  lifespan_acute = 0.5
)
model$set_parameters(
  acute_duration = 0.40,
  lifespan_acute = 1.4
)
model$update(1*400)
model$results_wide |> slice(360:370)
model$results_wide |>
  pivot_longer(S:Rf, values_to="Koalas", names_to="Compartment") |>
  ggplot(aes(x=Day, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(75), lty="dotted") +
  geom_vline(xintercept=365, lty="dashed")


## Observed population numbers:
tribble(
  ~Date, ~Outcome, ~Source, ~LCI, ~UCI,
  "2022-07-01", 300, "Chad", 281, 569,
  "2025-07-01", 107, "Lachlan", 95, 123,
) |>
  mutate(Date = as_date(Date)) |>
  mutate(Compartment = "Total Koalas (N)") |>
  identity() ->
  population

## Observed prevalence estimates:
tribble(
  ~Date, ~Positive, ~Total, ~Source,
  "2022-07-01", 2, 18, "Swab",
  "2022-11-01", 4, 17, "Swab",
  "2023-02-01", 5, 16, "Swab",
  "2023-07-01", 3, 13, "Swab",
  "2023-10-01", 4, 10, "Swab",
  "2021-07-01", 5, 11, "Scat",
  "2023-07-01", 10, 14, "Scat",
  "2025-07-01", 2, 16, "Scat",
) |>
  mutate(Date = as_date(Date), Outcome = Positive/Total*100) |>
  mutate(LCI = 100*qbeta(0.025, Positive+1, (Total-Positive)+1)) |>
  mutate(UCI = 100*qbeta(0.975, Positive+1, (Total-Positive)+1)) |>
  mutate(Compartment = "Observed Prevalence (%)") |>
  identity() ->
  prevalence


sensitivity <- KoalasV2$new()$parameters$sensitivity

pdf("calibration.pdf", width=7, height=5)

## Calibration for worst case scenario:
model <- KoalasV2$new(start_date = "2021-01-01")

prev <- 0.009
N <- 220
model$set_state(
  S = N * (1.0-prev),
  I = N * prev
)
model$set_parameters(
  acute_duration = 0.40,
  lifespan_acute = 1.4,
  beta = 2.69
)

model$update(365*5)
model$results_long |> filter(Date=="2022-07-01", Compartment=="Total") |> pull(Koalas)
population |> filter(Date=="2022-07-01")

model$results_long |>
  filter(Compartment %in% c("Infectious","Total")) |>
  mutate(Outcome = case_when(
    Compartment=="Total" ~ Koalas,
    .default = Percent * sensitivity
  ), Source = "Model") |>
  mutate(Compartment = case_when(
    Compartment=="Infectious" ~ "Observed Prevalence (%)",
    Compartment=="Total" ~ "Total Koalas (N)",
  )) |>
  ggplot(aes(x=Date, y=Outcome)) +
  geom_line() +
  facet_wrap(~Compartment, scales="free_y", nrow=2) +
  geom_point(mapping=aes(col=Source),
    data = bind_rows(prevalence, population),
    size = 2.5
  ) +
  geom_errorbar(mapping=aes(ymin=LCI, ymax=UCI, col=Source),
    data=bind_rows(prevalence, population),
    lty="dashed", width=75
  ) +
  ylim(c(0,NA)) +
  ylab(NULL) + xlab(NULL) +
  ggtitle(label="Worst case scenario")


## Calibration for best case scenario:
model <- KoalasV2$new(start_date="2021-01-01")
prev <- 0.009
N <- 61
model$set_state(
  S = N * (1.0-prev),
  I = N * prev
)
model$set_parameters(
  acute_duration = 0.40,
  lifespan_acute = 1.4,
  beta = 2.69
)
model$update(365*5)

model$results_long |> filter(Date=="2025-07-01", Compartment=="Total") |> pull(Koalas)
population |> filter(Date=="2025-07-01")

model$results_long |>
  filter(Compartment %in% c("Infectious","Total")) |>
  mutate(Outcome = case_when(
    Compartment=="Total" ~ Koalas,
    .default = Percent * sensitivity
  ), Source = "Model") |>
  mutate(Compartment = case_when(
    Compartment=="Infectious" ~ "Observed Prevalence (%)",
    Compartment=="Total" ~ "Total Koalas (N)",
  )) |>
  ggplot(aes(x=Date, y=Outcome)) +
  geom_line() +
  facet_wrap(~Compartment, scales="free_y", nrow=2) +
  geom_point(mapping=aes(col=Source),
    data = bind_rows(prevalence, population),
    size = 2.5
  ) +
  geom_errorbar(mapping=aes(ymin=LCI, ymax=UCI, col=Source),
    data=bind_rows(prevalence, population),
    lty="dashed", width=75
  ) +
  ylim(c(0,NA)) +
  ylab(NULL) + xlab(NULL) +
  ggtitle(label="Best case scenario")


## Calibration for hitting both (combined scenario):
model <- KoalasV2$new(start_date="2021-01-01")
prev <- 0.011
N <- 305
model$set_state(
  S = N * (1.0-prev),
  I = N * prev
)
model$set_parameters(
  acute_duration = 0.5,  # 0.4
  lifespan_acute = 0.5,  # 1.4
  lifespan_natural = 6,  # 6
  lifespan_chronic = 1.5,  # 4.8
  birthrate = 0.17,  # 0.38
  subclinical_recover_proportion = 0,  # 0.35
  beta = 2.22  # 2.64
)
model$update(365*5)

model$results_long |>
  inner_join(population |> select(Date, Outcome),
    join_by(Date)) |>
  filter(Compartment=="Total") |>
  mutate(Diff = Outcome-Koalas) |>
  pull(Diff)

model$results_long |>
  filter(Compartment %in% c("Infectious","Total")) |>
  mutate(Outcome = case_when(
    Compartment=="Total" ~ Koalas,
    .default = Percent * sensitivity
  ), Source = "Model") |>
  mutate(Compartment = case_when(
    Compartment=="Infectious" ~ "Observed Prevalence (%)",
    Compartment=="Total" ~ "Total Koalas (N)",
  )) |>
  ggplot(aes(x=Date, y=Outcome)) +
  geom_line() +
  facet_wrap(~Compartment, scales="free_y", nrow=2) +
  geom_point(mapping=aes(col=Source),
    data = bind_rows(prevalence, population),
    size = 2.5
  ) +
  geom_errorbar(mapping=aes(ymin=LCI, ymax=UCI, col=Source),
    data=bind_rows(prevalence, population),
    lty="dashed", width=75
  ) +
  ylim(c(0,NA)) +
  ylab(NULL) + xlab(NULL) +
  ggtitle(label="Combined scenario")

dev.off()
