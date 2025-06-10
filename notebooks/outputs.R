## Outputs for report

library("koalas")
library("tidyverse")
theme_set(theme_light())

model <- KoalasV2$new()

model$set_state(S = 290, I = 10, N=100)

model$set_parameters(
  passive_intervention_rate = 0
)
model$update(1*365)

model$set_parameters(
  passive_intervention_rate = 0.01,
  sensitivity = 1,
  specificity = 1
)
model$update(1*365)

for(y in seq_len(9)){
  model$update(182)
  model$active_intervention(prop=1, cull_positive=0, cull_acute=0.2, cull_chronic=0.3)
  model$update(365-182)
}

model$results_wide |> tail() |> View()

model$results_long |>
  mutate(Date = as.Date("2024-01-01") + Year*365) |>
  ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(11,30), lty="dotted")



## TODO: calibrate acute mortality to 75% survival

## Destinations after treatment:
# N:  sp are missed, 0% -> Z, of remaining: sp go to N, of remaining: 90% -> R, 10% -> N
# I:  1-se are missed, 0% -> Z, of remaining: 1-se go to I, of remaining: 90% -> R, 10% -> N
# A:  0% are missed, 20% -> Z, of remaining: 1-se go to I, of remaining: 75% -> R, 25% -> N
# C:  0% are missed, 30% -> Z, of remaining: 1-se go to I, of remaining: 60% -> R, 40% -> N

# S:  1-sp are treated as I, remaining %->V, ->S

## Scenarios - goal: keeping prevalence under 5% vs 10% vs 20%, for 2 vs 4 events per year (single)
# 0. No intervention, just passive at 3%
# Start with e.g. 85% sampling for 2 years, then reduce to e.g. 25% sampling
# 1. Soft approach: minimum sampling coverage to go down to 5% then maintain it (e.g. start with 85%)
# 2. Harder approach: increase culling of chronic from 30% to 50%, 75%, 100%
# 3. Even harder approach: increase culling of acute+chronic from 20/30% to 50%, 75%, 100%
# 4. Even harder approach: increase culling of infected+acute+chronic from 0/20/30% to 50%, 75%, 100%
# Present:
# 1. prevalence % (I+If+Af+Cf)
# 2. total number
# 3. cumulative treatments (not including %culled, i.e. not including failed treatments)


## passive number: 3%

## Start: 11% in mid 2022, 30% in mid 2023, ?? in mid 2025
## Interventions:  Q4 2025, Oct+April, plus Jan+July

prev <- 0.11
N <- 300
model$set_state(
  S = N * (1.0-prev),
  V = 0.0,
  I = N * prev * 0.6,
  N = 0.0,
  R = 0.0,
  Af = N * prev * 0.3,
  Cf = N * prev * 0.1,
  Sf = 0.0,
  Vf = 0.0,
  If = 0.0,
  Nf = 0.0,
  Rf = 0.0
)

model$set_parameters(
  beta = 2.1,
  #lifespan_acute = 0.5,  # 25% die before moving to chronic within 1 year
  lifespan_chronic = 4.8,  # 80% of natural life
  lifespan_natural = 6,   #
  passive_intervention_rate = 0,
  treatment_dest_R = 1,
  treatment_dest_N = 0,
  treatment_dest_IAC = 0,
  treatment_dest_remove = 0
)

model$update(1*365)

model$set_parameters(
  #beta = 2.1,
  #lifespan_acute = 0.5,
  #lifespan_chronic = 1,
  passive_intervention_rate = 0.01
)

for(y in seq_len(9)){
  model$update(182)
  model$active_intervention(prop=1, cull_positive=0.9, cull_acute=0.9, cull_chronic=0.9)
  model$update(365-182)
}

model$results_wide

model$results_long |>
  mutate(Date = as.Date("2024-01-01") + Year*365) |>
  ggplot(aes(x=Date, y=Percent, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(11,30), lty="dotted")

model$results_long |>
  mutate(Date = as.Date("2024-01-01") + Year*365) |>
  ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
  geom_line()



model$set_parameters(
  passive_intervention_rate = 0.01
)


ggplot(model$results_long, aes(x=Year, y=Koalas, col=Compartment)) +
  geom_line()

model$results_long |>
  mutate(Date = as.Date("2022-01-01") + Year*365) |>
  ggplot(aes(x=Date, y=Percent, col=Compartment)) +
  geom_line() +

