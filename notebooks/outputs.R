## Outputs for report

library("koalas")
library("tidyverse")
theme_set(theme_light())

## Outputs
# 1. Baseline from July 2022 to October 2035, plus dotted line for 2022+2023, 3 dashed lines
# 2. Sensitivity graphs for 2 years
# 3. Sensitivity graphs for 8 years
# 4. 3-row graphs for: initial=4 x ~45%, later=2 x ~40%, 4 dashed lines

subfolder <- "reports/partialvacc/"
if(!dir.exists(subfolder)) dir.create(subfolder)

## NOTE: when lifespan_natural < Inf vaccination decreases population size
## R, Rf, V and Vf animals die faster than S and Sf

## 0. Current prevalence is expected to be around 65%
model <- KoalasV2$new()
model$set_state(S=300, V=0)
model$set_parameters(vaccine_efficacy=0, vaccine_booster=0, lifespan_natural=6, passive_intervention_rate=0, birthrate=0)
sum(unlist(model$state[1:12]))
model$run(10, frequency = 0)
sum(unlist(model$state[1:12]))

model$parameters
model$set_parameters(vaccine_efficacy=0, vaccine_booster=0, lifespan_natural=Inf)
model$burnin()
model$results_long |> filter(Date=="2025-07-01")
model$results_long |>
  filter(Compartment == "Infectious") |>
  ggplot(aes(x=Date, y=Percent)) +
  geom_line() +
  geom_vline(xintercept=as.Date(c("2022-07-01","2023-07-01","2025-07-01")), lty="dashed") +
  geom_hline(yintercept=c(10,30,65), lty="dotted") +
  ylim(0,100) +
  ylab("Prevalence")
# We also hit the targets of 10% and 30% as expected
ggsave(file.path(subfolder, "Figure 0.pdf"), height=5, width=6)


## 1. Baseline
baseline <- model$clone(deep=TRUE)
baseline$run(10, frequency = 0)
baseline$autoplot() +
  geom_vline(xintercept=as.Date(c("2022-07-01","2023-07-01","2025-07-01")), lty="dashed")
ggsave(file.path(subfolder, "Figure 1.pdf"), height=6, width=6)

# Or maybe:
baseline$autoplot(show_treatments=FALSE) +
  geom_vline(xintercept=as.Date(c("2022-07-01","2023-07-01","2025-07-01")), lty="dashed")


## 2. Initial phase
phase1 <- assess_interventions(
  model,
  years = 2,
  frequency = 1:4,
  prop_active = seq(0, 1, by=0.01),
  prop_targeted = 0,
  cl = 10L  # Note: remove this argument if you see error messages
)
phase1 |>
  filter(Prevalence <= 5) |>
  group_by(Frequency) |>
  arrange(PropActive) |>
  slice(1L) ->
  optimums
optimums
# We end up with 1=>100%, 2=76%, 3=58%, 4=47%

phase1 |>
  pivot_longer("Prevalence":"Koalas", names_to="Metric", values_to="Value") |>
  ggplot(aes(x=PropActive, y=Value, col=Frequency)) +
  geom_line() +
  geom_hline(data=tibble(Metric="Prevalence", ll=5), mapping=aes(yintercept=ll), lty="dashed") +
  geom_vline(data=optimums, mapping=aes(xintercept=PropActive, col=Frequency), lty="dashed") +
  facet_wrap(~Metric, ncol=1, scales="free_y") +
  theme(legend.position="right") +
  ylab(NULL) + xlab("Active Sampling Proportion") +
  theme(
    strip.background = element_rect(fill="grey95"),
    strip.text=element_text(color="black")
  )
ggsave(file.path(subfolder, "Figure 2.pdf"), height=4, width=6)


## 2b. Early phase with targeted interventions
phase1 <- assess_interventions(
  model,
  years = 2,
  frequency = 1:4,
  prop_active = seq(0, 1, by=0.01),
  prop_targeted = c(0,0.25,0.5,0.75),
  cl = 10L  # Note: remove this argument if you see error messages
)
phase1 |>
  filter(Prevalence <= 5) |>
  group_by(Frequency, PropTargeted) |>
  arrange(PropActive) |>
  slice(1L) ->
  optimums
optimums

phase1 |>
  ggplot(aes(x=PropActive, y=Prevalence, col=Frequency)) +
  geom_line() +
  geom_hline(yintercept=5, lty="dashed") +
  geom_vline(data=optimums, mapping=aes(xintercept=PropActive, col=Frequency), lty="dashed") +
  facet_wrap(~str_c("Targeted Sampling: ", PropTargeted), ncol=1, scales="free_y") +
  theme(legend.position="right") +
  ylab(NULL) + xlab("Active Sampling Proportion") +
  theme(
    strip.background = element_rect(fill="grey95"),
    strip.text=element_text(color="black")
  ) +
  ylim(0,100)
ggsave(file.path(subfolder, "Figure 2b.pdf"), height=8, width=6)

## 2c. Early phase with only targeted interventions
phase1 <- assess_interventions(
  model,
  years = 2,
  frequency = 1:4,
  prop_active = c(0, 0.25, 0.5, 0.75),
  prop_targeted = seq(0, 1, by=0.01),
  cl = 10L  # Note: remove this argument if you see error messages
)
phase1 |>
  filter(Prevalence <= 5) |>
  group_by(Frequency, PropActive) |>
  arrange(PropTargeted) |>
  slice(1L) ->
  optimums
optimums

phase1 |>
  ggplot(aes(x=PropTargeted, y=Prevalence, col=Frequency)) +
  geom_line() +
  geom_hline(yintercept=5, lty="dashed") +
  geom_vline(data=optimums, mapping=aes(xintercept=PropTargeted, col=Frequency), lty="dashed") +
  facet_wrap(~str_c("Active Sampling: ", PropActive), ncol=1, scales="free_y") +
  theme(legend.position="right") +
  ylab(NULL) + xlab("Targeted Sampling Proportion") +
  theme(
    strip.background = element_rect(fill="grey95"),
    strip.text=element_text(color="black")
  ) +
  ylim(0,100)
ggsave(file.path(subfolder, "Figure 2c.pdf"), height=8, width=6)


## 3. Later phase
model$run(2, frequency = 4, prop_active = 0.47, prop_targeted = 0.0)

phase2 <- assess_interventions(
  model,
  years = 8,
  frequency = 1:4,
  prop_active = seq(0, 1, by=0.01),
  prop_targeted = 0.0,
  cl = 10L  # Note: remove this argument if you see error messages
)
phase2 |>
  filter(Prevalence <= 5) |>
  group_by(Frequency) |>
  arrange(PropActive) |>
  slice(1L) ->
  optimums
optimums
# We end up with 1=>75%, 2=45%, 3=32%, 4=25%

phase2 |>
  pivot_longer("Prevalence":"Koalas", names_to="Metric", values_to="Value") |>
  ggplot(aes(x=PropActive, y=Value, col=Frequency)) +
  geom_line() +
  geom_hline(data=tibble(Metric="Prevalence", ll=5), mapping=aes(yintercept=ll), lty="dashed") +
  geom_vline(data=optimums, mapping=aes(xintercept=PropActive, col=Frequency), lty="dashed") +
  facet_wrap(~Metric, ncol=1, scales="free_y") +
  theme(legend.position="right") +
  ylab(NULL) + xlab("Active Sampling Proportion") +
  theme(
    strip.background = element_rect(fill="grey95"),
    strip.text=element_text(color="black")
  )
ggsave(file.path(subfolder, "Figure 3.pdf"), height=4, width=6)


## 3b. With targeted surveillance

phase2 <- assess_interventions(
  model,
  years = 8,
  frequency = 1:4,
  prop_active = seq(0, 1, by=0.01),
  prop_targeted = c(0,0.25,0.5,0.75),
  cl = 10L  # Note: remove this argument if you see error messages
)
phase2 |>
  filter(Prevalence <= 5) |>
  group_by(Frequency, PropTargeted) |>
  arrange(PropActive) |>
  slice(1L) ->
  optimums
optimums

phase2 |>
  ggplot(aes(x=PropActive, y=Prevalence, col=Frequency)) +
  geom_line() +
  geom_hline(yintercept=5, lty="dashed") +
  geom_vline(data=optimums, mapping=aes(xintercept=PropActive, col=Frequency), lty="dashed") +
  facet_wrap(~str_c("Targeted Sampling: ", PropTargeted), ncol=1, scales="free_y") +
  theme(legend.position="right") +
  ylab(NULL) + xlab("Active Sampling Proportion") +
  theme(
    strip.background = element_rect(fill="grey95"),
    strip.text=element_text(color="black")
  ) +
  ylim(0,100)
ggsave(file.path(subfolder, "Figure 3b.pdf"), height=8, width=6)



## 4. Final scenario
model$run(8, frequency = 2, prop_active = 0.45, prop_targeted = 0)
model$autoplot() +
  geom_vline(xintercept=model$run_dates, lty="dashed")
ggsave(file.path(subfolder, "Figure 4.pdf"), height=6, width=6)
