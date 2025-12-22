
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



## Calibration of beta and burnin phase starting at 5% and ending at 10% with 300 koalas then 30% 1 year later
model <- KoalasV2$new()
prev <- 0.05
N <- 270
model$set_state(
  S = N * (1.0-prev),
  I = N * prev
)
model$set_parameters(
  subclinical_recover_proportion = 0.35,
  beta = 2.35
)
dd <- 175
model$update(dd+400)

model$results_long |>
  mutate(Day = as.numeric(Date - min(Date), units="days")) |>
  ggplot(aes(x=Day, y=Percent, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(10,30), lty="dotted") +
  geom_vline(xintercept=c(dd,dd+365), lty="dashed") +
  labs(title=model$parameters$subclinical_recover_proportion)

model$results_long |>
  filter(Compartment=="Total") |>
  mutate(Day = as.numeric(Date - min(Date), units="days")) |>
  ggplot(aes(x=Day, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(300), lty="dotted") +
  geom_vline(xintercept=c(dd), lty="dashed") +
  labs(title=model$parameters$subclinical_recover_proportion)


