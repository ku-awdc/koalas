## Outputs for report

library("koalas")
library("tidyverse")
theme_set(theme_light())
library("pbapply")

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


expand_grid(
  Frequency = c(2,4),
  InitialProp = seq(0.5,0.9,by=0.1),
  LaterProp = seq(0.1,0.5,by=0.1),
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

pdf("output.pdf", height=10, width=8)
scenarios |>
#  slice_sample(n=10) |>
  rowwise() |>
  group_split() |>
  pblapply(\(x){

    mm <- KoalasV2$new()
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

