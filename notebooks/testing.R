library("koalas")
library("tidyverse")
library("pbapply")
theme_set(theme_light())

expand_grid(
  Scenario = letters[1:6] |> fct(levels = letters[1:6] |> matrix(nrow=2) |> t() |> as.character()),
  Capacity = c("Unlimited","Limited")
) |>
  mutate(
    carrying_capacity = if_else(Capacity=="Unlimited", Inf, 1000),
    S = 800,
    I = 200,
    beta = 3.0,
    sigma = 0.6,
    acute_recovery_prob = 0.29,
    birth_rate = case_when(
      Scenario=="b" ~ 0.15,
      TRUE ~ 0.38
    ),
    mortality_natural = case_when(
      Scenario=="b" ~ 0.3,
      TRUE ~ 0.2
    ),
    mortality_disease = case_when(
      Scenario%in%c("c","e","f") ~ 0.5,
      Scenario%in%c("d") ~ 0.9,
      Scenario%in%c("a","b") ~ 0,
    ),
    relative_fecundity = case_when(
      Scenario%in%c("c","e","f") ~ 0.5,
      Scenario%in%c("d") ~ 0.1,
      Scenario%in%c("a","b") ~ 1.0,
    ),
    recovery = case_when(
      Scenario%in%c("e","f") ~ 0.1,
      TRUE ~ 0
    ),
    waning_natural = case_when(
      Scenario%in%c("f") ~ 0.5,
      TRUE ~ 0
    ),
    waning_vaccine = 0
  ) ->
  parameters

if(FALSE){
  parameters |>
    filter(Scenario=="a", Capacity=="Limited") ->
    x
  model <- KoalasR6$new()
  for(nn in names(x)){
    if(nn %in% c("Scenario","Capacity")) next
    model[[nn]] <- x[[nn]]
  }
  model$run(20, 0.01) |>
    select(Time:R, N) |>
    pivot_longer(-Time, names_to="Compartment", values_to="Number") |>
    bind_cols(x |> select(Scenario, Capacity))
}

parameters |>
  rowwise() |>
  group_split() |>
  pblapply(\(x){
    model <- KoalasR6$new()
    for(nn in names(x)){
      if(nn %in% c("Scenario","Capacity")) next
      model[[nn]] <- x[[nn]]
    }
    model$run(20, 0.01) |>
      select(Time:R, N) |>
      pivot_longer(-Time, names_to="Compartment", values_to="Number") |>
      bind_cols(x |> select(Scenario, Capacity))
  }, cl=6) |>
  bind_rows() ->
  res

## Without CC
res |>
  filter(Capacity=="Unlimited") |>
  ggplot(aes(x=Time, y=Number, col=Compartment)) +
  geom_line() +
  ylim(0,1500) +
  facet_wrap(~Scenario) +
  scale_color_manual(values=c(S="forestgreen", I="red", D="purple", R="blue", N="grey50")) +
  ylab("Number of koalas") + xlab("Time (years)")
ggsave("notebooks/fig2.pdf", width=12, height=7)

## With CC
res |>
  filter(Capacity=="Limited") |>
  ggplot(aes(x=Time, y=Number, col=Compartment)) +
  geom_line() +
  ylim(0,1500) +
  facet_wrap(~Scenario) +
  scale_color_manual(values=c(S="forestgreen", I="red", D="purple", R="blue", N="grey50")) +
  ylab("Number of koalas") + xlab("Time (years)")
ggsave("notebooks/figS3.pdf", width=12, height=7)

