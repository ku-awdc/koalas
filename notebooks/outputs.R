## Outputs for report

library("koalas")
library("tidyverse")
theme_set(theme_light())

mm <- KoalasV2$new()
mm$burnin()
mm$date

mm$results_long |>
  ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(300), lty="dotted")

mm$results_long |>
  ggplot(aes(x=Date, y=Percent, col=Compartment)) +
  geom_line() +
  geom_hline(yintercept=c(10,30), lty="dotted") +
  geom_vline(xintercept=as.Date("2022-07-01")+365*c(0,1,3), lty="dashed")

mm$run(years = 2, sampling_frequency = 2, proportion = 0.75)

mm$results_long |>
  ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
  geom_line()
mm$results_long |>
  ggplot(aes(x=Date, y=Percent, col=Compartment)) +
  geom_line()

mm$run(years = 8, sampling_frequency = 2, proportion = 0.5)
mm$results_long |>
  ggplot(aes(x=Date, y=Koalas, col=Compartment)) +
  geom_line()
mm$results_long |>
  ggplot(aes(x=Date, y=Percent, col=Compartment)) +
  geom_line()

mm$treatments |>
  ggplot(aes(x=Date, y=Cumulative, col=Type)) +
  geom_line()
