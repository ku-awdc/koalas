
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
  acute_duration = 0.40,
  lifespan_acute = 4.0
)
model$update(1*400)
model$results_wide |> slice(360:370)
model$results_wide |>
  pivot_longer(S:Rf, values_to="Koalas", names_to="Compartment") |>
  ggplot(aes(x=Day, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(75), lty="dotted") +
  geom_vline(xintercept=365, lty="dashed")


## Calibration of beta so that we go from 10% to 30% prevalence in 1 year
model <- KoalasV2$new()
prev <- 0.10
N <- 300
model$set_state(
  S = N * (1.0-prev),
  I = N * prev * 0.6,
  Af = N * prev * 0.3,
  Cf = N * prev * 0.1
)
model$set_parameters(
  beta = 2.2
)
model$update(400)
model$results_long |>
  ggplot(aes(x=Year, y=Percent, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(10,30), lty="dotted") +
  geom_vline(xintercept=1, lty="dashed")


## Calibration of burnin phase starting at 5% and ending at 10% with 300 koalas
model <- KoalasV2$new()
prev <- 0.05
N <- 272
model$set_state(
  S = N * (1.0-prev),
  I = N * prev
)
model$set_parameters(
  beta = 2.2
)
model$update(400)
dd <- 197
model$results_long |>
  ggplot(aes(x=Year, y=Percent, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=10, lty="dashed") +
  geom_vline(xintercept=dd/365, lty="dotted")
model$results_long |>
  ggplot(aes(x=Year, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=300, lty="dashed") +
  geom_vline(xintercept=dd/365, lty="dotted")

