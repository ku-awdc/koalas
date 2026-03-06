## Outputs for report

library("koalas")
library("tidyverse")
library("lubridate")
theme_set(theme_light())


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
  filter(Source=="Swab") |>
  mutate(Date = as_date(Date), Outcome = (1/0.95)*Positive/Total*100) |>
  mutate(LCI = 100*(1/0.95)*qbeta(0.025, Positive+1, (Total-Positive)+1)) |>
  mutate(UCI = 100*(1/0.95)*qbeta(0.975, Positive+1, (Total-Positive)+1)) |>
  mutate(Compartment = "Observed Prevalence (%)", Percent=Outcome) |>
  identity() ->
  prevalence


## Outputs
# 1. Baseline from July 2022 to October 2035, plus dotted line for 2022+2023, 3 dashed lines
# 2. Sensitivity graphs for 2 years
# 3. Sensitivity graphs for 8 years
# 4. 3-row graphs for: initial=4 x ~45%, later=2 x ~40%, 4 dashed lines

## Note: we were previously over-estimating mortality by 3x for V/N/I/R/A compartments
scenario <- "best"
scenario <- "worst"
#scenario <- "combined"

with_obs <- TRUE
with_ci <- TRUE

if(scenario=="worst"){
  pop <- population |> filter(Source=="Chad")
}else if(scenario=="best"){
  pop <- population |> filter(Source=="Lachlan")
}else{
  pop <- population
}

prevalence <- prevalence |> mutate(Subplot = "Prevalence of Infection (%)")
pop <- pop |> mutate(Subplot = "Number of Koalas")

months <- as_date("2022-07-01") + months(0:(12*4))
quarters <- as_date("2022-07-01") + months(0:(12*4/3)*3)
years <- as_date("2022-01-01") + years(0:15)

nn <- str_c("secnario_", scenario)
subfolder <- file.path("reports",nn)
if(!dir.exists(subfolder)) dir.create(subfolder)

## 0. Current prevalence is expected to be around 65%
model <- KoalasV2$new()
#model$parameters

## Remove vaccination effect:
if(nn=="novacc"){
  model$set_parameters(vaccine_efficacy=0, vaccine_booster=0)
}

model$burnin(scenario=scenario)
model$results_long |> filter(Date=="2025-07-01")
model$results_long |>
  filter(Compartment == "Infectious") |>
  ggplot(aes(x=Date, y=Percent)) +
  geom_line() +
  geom_vline(xintercept=as.Date(c("2026-06-01")), lty="dashed") +
#  geom_hline(yintercept=c(10,30,65), lty="dotted") +
  ylim(0,100) +
  ylab("Prevalence") + xlab(NULL) +
  scale_x_date(date_labels="%Y", minor_breaks=quarters) +
  geom_point(
    data = prevalence,
    size = 2.5
  ) +
  geom_errorbar(mapping=aes(ymin=LCI, ymax=UCI),
    data=prevalence,
    lty="dashed", width=75
  )
ggsave(file.path(subfolder, "Figure 0.pdf"), height=5, width=6)


## 1. Baseline
baseline <- model$clone(deep=TRUE)
baseline$run(10, frequency = 0)
pt <- baseline$autoplot() +
  # geom_vline(xintercept=as.Date(c("2022-07-01","2023-07-01","2025-07-01")), lty="dashed") +
  scale_x_date(minor_breaks=years)

if(with_obs){
  pt <- pt +
    geom_point(
      mapping = aes(x=Date, y=Outcome),
      data = bind_rows(prevalence, pop),
      size = 1.5
    )
}
if(with_ci){
  pt <- pt +
    geom_errorbar(
      mapping = aes(x=Date, y=Outcome, ymin=LCI, ymax=UCI),
      data = bind_rows(prevalence, pop),
      lty="dotted", width=0
    )
}
pt #+ geom_vline(xintercept=as.Date(c("2026-06-01")), lty="dashed")
ggsave(file.path(subfolder, "Figure 1.pdf"), height=6, width=6)

# Or maybe:
baseline$autoplot(show_treatments=FALSE) +
  geom_vline(xintercept=as.Date(c("2022-07-01","2023-07-01","2025-07-01")), lty="dashed")

# Extract output:
baseline$results_wide |>
  select("Date":"Rf") |>
  mutate(Total = rowSums(across(-c("Date", "Day")))) |>
  filter(row_number()<=2L | .data$Total > 0) |> # Remove rows with NA prevalence, but keep first 2 rows to make sure the plot is created
  mutate(Sum = .data$I+.data$If+.data$Af+.data$Cf, Prevalence = .data$Sum/.data$Total * 100) |>
  select("Date", "Population"="Total", "Prevalence") |>
  filter(day(Date)==1) |>
  left_join(
    baseline$treatments |> pivot_wider(names_from="Type", values_from="Cumulative"),
    join_by(Date)
  ) |>
  identity() ->
  baseline_output


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
pt <- model$autoplot() +
  geom_vline(xintercept=model$run_dates, lty="dashed") +
  scale_x_date(minor_breaks=years)
if(with_obs){
  pt <- pt +
    geom_point(
      mapping = aes(x=Date, y=Outcome),
      data = bind_rows(prevalence, pop),
      size = 1.5
    )
}
if(with_ci){
  pt <- pt +
    geom_errorbar(
      mapping = aes(x=Date, y=Outcome, ymin=LCI, ymax=UCI),
      data = bind_rows(prevalence, pop),
      lty="dotted", width=0
    )
}
#pt + geom_vline(xintercept=as.Date(c("2026-06-01")), lty="dashed")
ggsave(file.path(subfolder, "Figure 4.pdf"), height=6, width=6)

model$results_wide |>
  select("Date":"Rf") |>
  mutate(Total = rowSums(across(-c("Date", "Day")))) |>
  filter(row_number()<=2L | .data$Total > 0) |> # Remove rows with NA prevalence, but keep first 2 rows to make sure the plot is created
  mutate(Sum = .data$I+.data$If+.data$Af+.data$Cf, Prevalence = .data$Sum/.data$Total * 100) |>
  select("Date", "Population"="Total", "Prevalence") |>
  filter(day(Date)==1) |>
  left_join(
    model$treatments |> pivot_wider(names_from="Type", values_from="Cumulative"),
    join_by(Date)
  ) |>
  identity() ->
  final_output

writexl::write_xlsx(
  list(baseline=baseline_output, interventions=final_output),
  file.path(subfolder, str_c("output_", scenario, ".xlsx"))
)

ff <- str_c("Figure ", c("0","1","2","2b","2c","3","3b","4"), ".pdf")
(cmd <- str_c("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=", file.path(subfolder,str_c(nn,".pdf")), " ", str_c("'", file.path(subfolder,ff), "'", collapse=" ")))
system(cmd)

# gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=out.pdf in1.pdf in2.pdf

