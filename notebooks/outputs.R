## Outputs for report

library("koalas")
library("tidyverse")
theme_set(theme_light())
library("pbapply")

## Outputs
# 1. Baseline from July 2022 to October 2035, plus dotted line for 2022+2023, 3 dashed lines
# 2. Sensitivity graphs for 2 years
# 3. Sensitivity graphs for 8 years
# 4. 3-row graphs for: initial=4 x ~45%, later=2 x ~40%, 4 dashed lines

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


# Optimisation

expand_grid(
  Frequency = c(2,4),
  InitialProp = seq(0.5,0.9,by=0.1),
  LaterProp = 0.25, #seq(0.1,0.5,by=0.1),
  CullPositive = c(0, 0.5, 1),
  CullAcute = c(0.2, 0.5, 1),
  CullChronic = c(0.3, 0.5, 1)
) |>
  bind_rows(
    tibble(Frequency = 0, InitialProp = 0, LaterProp = 0, CullPositive = 0, CullAcute = 0, CullChronic = 0)
  ) |>
  filter(CullChronic >= CullAcute, CullAcute >= CullPositive) ->
  scenarios

# x <- scenarios |> slice_sample(n=1)


findprop1 <- function(x, give_model=FALSE){

  mm <- KoalasV2$new()
  mm$burnin()
  mm$run(years = 2, sampling_frequency = x$Frequency, proportion = x$InitialProp, cull_positive = x$CullPositive, cull_acute = x$CullAcute, cull_chronic = x$CullChronic)

  if(give_model) return(mm)

  mm$results_wide |>
    slice(n()) |>
    select(Date:Rf) |>
    mutate(Total = rowSums(across(-c(Date, Day))), Sum=Total) |>
    mutate(Healthy=S+V+N+R, Infectious=I+If+Af+Cf, Diseased=Af+Cf, Infertile=Sf+Vf+Nf+Rf+If+Diseased, Immune=V+Vf+R+Rf+N+Nf) |>
    select(Date, Total, Infectious) ->
    final_prev

  x |>
    bind_cols(final_prev) |>
    mutate(Prevalence = Infectious / Total)
}

findprop2 <- function(x, give_model=FALSE){

  mm <- KoalasV2$new()
  mm$set_parameters(
    subclinical_recover_proportion = x$RecoverProp,
    beta = case_match(x$RecoverProp,
      0.05 ~ 2.25,
      0.35 ~ 2.6
    )
  )
  mm$burnin()
  mm$run(years = 2, sampling_frequency = 2, proportion = 0.75, cull_positive = x$CullPositive, cull_acute = x$CullAcute, cull_chronic = x$CullChronic)

  mm$run(years = 8, sampling_frequency = x$Frequency, proportion = x$LaterProp, cull_positive = x$CullPositive, cull_acute = x$CullAcute, cull_chronic = x$CullChronic)

  if(give_model) return(mm)

  mm$results_wide |>
    slice(n()) |>
    select(Date:Rf) |>
    mutate(Total = rowSums(across(-c(Date, Day))), Sum=Total) |>
    mutate(Healthy=S+V+N+R, Infectious=I+If+Af+Cf, Diseased=Af+Cf, Infertile=Sf+Vf+Nf+Rf+If+Diseased, Immune=V+Vf+R+Rf+N+Nf) |>
    select(Date, Total, Infectious) ->
    final_prev

  x |>
    bind_cols(final_prev) |>
    mutate(Prevalence = Infectious / Total)
}

expand_grid(
  Frequency = c(1,2,3,4),
  RecoverProp = 0.35, #c(0.05, 0.35),
  InitialProp = seq(0, 1, by=0.1),
  LaterProp = 0, #seq(0.0,0.75,by=0.05),
  CullPositive = 0, #
  CullAcute = 0.2, #c(0.2, 0.5, 1),
  CullChronic = 0.3 #c(0.3, 0.5, 1)
) |>
  rowwise() |>
  group_split() |>
  pblapply(findprop1, cl=10L) |>
  bind_rows() ->
  res

ggplot(res, aes(x=InitialProp, y=Prevalence*100, col=factor(Frequency))) +
  geom_line() +
  geom_hline(yintercept=5, lty="dashed") +
  ylim(c(0,NA_real_))

ggplot(res, aes(x=InitialProp, y=Total, col=factor(Frequency))) +
  geom_line() +
  ylim(c(0,NA_real_))


expand_grid(
  Frequency = c(1,2,3,4),
  RecoverProp = 0.35, #c(0.05, 0.35),
  LaterProp = seq(0, 1, by=0.1),
  CullPositive = 0, #
  CullAcute = 0.2, #c(0.2, 0.5, 1),
  CullChronic = 0.3 #c(0.3, 0.5, 1)
) |>
  rowwise() |>
  group_split() |>
  pblapply(findprop2, cl=10L) |>
  bind_rows() ->
  res

ggplot(res, aes(x=LaterProp, y=Prevalence*100, col=factor(Frequency))) +
  geom_line() +
  geom_hline(yintercept=5, lty="dashed") +
  ylim(c(0,NA_real_))

ggplot(res, aes(x=LaterProp, y=Total, col=factor(Frequency))) +
  geom_line() +
  ylim(c(0,NA_real_))


expand_grid(
  Frequency = 4, #c(1,2,3,4),
  RecoverProp = 0.35, #c(0.05, 0.35),
  InitialProp = 0,
  LaterProp = 0.25,
  CullPositive = 0, #
  CullAcute = 0.2, #c(0.2, 0.5, 1),
  CullChronic = 0.3 #c(0.3, 0.5, 1)
) -> x

model <- findprop1(x, TRUE)
model$autoplot()

model$results_long |>
  filter(Compartment=="Infectious") |>
  ggplot(aes(x=Date, y=Percent)) +
  geom_line() +
  geom_vline(xintercept=as.Date(c("2022-07-01","2023-07-01","2025-10-01", "2027-10-01")), lty="dashed", col="grey") +
  geom_hline(yintercept=5, lty="dotted")

## Stacked bar chart, Diseased, Infected (I+If), Other
model$results_long |>
  filter(Compartment%in%c("Diseased","Infectious","Total")) |>
  ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_vline(xintercept=as.Date(c("2025-10-01", "2027-10-01")), lty="dashed") +
  ylim(c(0,NA_real_))

model$treatments |>
  filter(Type=="Treated") |>
  ggplot(aes(x=Date, y=Cumulative)) +
  geom_line()

expand_grid(
  Frequency = c(2,3,4,6,8),
  RecoverProp = c(0.05, 0.35),
  InitialProp = seq(0.2,1.0,by=0.1),
  LaterProp = 0.25, #seq(0.1,0.5,by=0.1),
  CullPositive = 0, #c(0, 0.5, 1),
  CullAcute = 0.2, #c(0.2, 0.5, 1),
  CullChronic = 0.3 #c(0.3, 0.5, 1)
) |>
  rowwise() |>
  group_split() |>
  pblapply(findprop1) |>
  bind_rows() ->
  res

ggplot(res, aes(x=InitialProp, y=FinalPrev, col=factor(Frequency))) +
  geom_line() +
  geom_hline(yintercept=5, lty="dashed") +
  ylim(c(0,NA_real_)) +
  facet_wrap(~RecoverProp, ncol=1)

pdf("output.pdf", height=10, width=8)
scenarios |>
#  slice_sample(n=10) |>
  rowwise() |>
  group_split() |>
  pblapply(\(x){

    mm <- KoalasV2$new()
    model$set_parameters(
      subclinical_recover_proportion = 0.35,
      beta = 2.6
    )
    mm$burnin()
    mm$run(years = 2, sampling_frequency = x$Frequency, proportion = x$InitialProp, cull_positive = x$CullPositive, cull_acute = x$CullAcute, cull_chronic = x$CullChronic)
    mm$run(years = 8, sampling_frequency = x$Frequency, proportion = x$LaterProp, cull_positive = x$CullPositive, cull_acute = x$CullAcute, cull_chronic = x$CullChronic)

    p1 <- mm$results_long |>
      filter(Compartment %in% c("Total","Infectious")) |>
      ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
      geom_line()
    p2 <- mm$results_long |>
      filter(Compartment %in% "Infectious") |>
      ggplot(aes(x=Date, y=Percent, col=Compartment)) +
      geom_line()
    p3 <- mm$treatments |>
      ggplot(aes(x=Date, y=Cumulative, col=Type)) +
      geom_line()

    tt <- str_c("F: ", x$Frequency, ", IP: ", x$InitialProp, ", LP: ",  x$LaterProp, ", CP: ", x$CullPositive, ", CA: ", x$CullAcute, ", CC: ", x$CullChronic)
    gridExtra::grid.arrange(p1, p2, p3, top=tt)

    invisible(TRUE)
  })
dev.off()

