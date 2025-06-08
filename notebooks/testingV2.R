library("tidyverse")
library("koalas")

## Compare to Grogan et al
mm <- KoalasV2$new(3)
mm$parameters$acute_duration
mm$set_parameters(acute_duration = 10)
mm$parameters$acute_duration
mm$state$R
mm$set_state(R=10, Day=1)
mm$state$R
mm$state$Day

mm$N
mm$state
mm$parameters
mm$results_wide
mm$update(50*365)


mm$results_wide |>
  select(Year,Day,S,I,R,Af,Cf) |>
  mutate(Total = S+I+R+Af+Cf) |>
  pivot_longer(S:Total) |>
  ggplot(aes(x=Year + Day/365, y=value, col=name)) +
  geom_line() +
  labs(caption=str_c("beta: ", pars[1], ",  birthrate: ", pars[2]), x="Time (years)", y="Number") +
  scale_x_continuous(breaks=0:10) +
  geom_vline(xintercept=1, lty="dashed")

mm$set_state(V=10)

mm <- KoalasV2$new(1)



pp <- mm$.__enclos_env__$private$.obj$parameters
pp["beta"] <- pars[1]
pp["birthrate"] <- pars[2]
pp["passive_intervention_rate"] <- 0.5
mm$.__enclos_env__$private$.obj$parameters <- pp




ff <- function(pars=c(3.0, 0.38)){
  mm <- KoalasV2$new()
  pp <- mm$.__enclos_env__$private$.obj$pars_natural
  pp["beta"] <- pars[1]
  pp["birthrate"] <- pars[2]
  mm$.__enclos_env__$private$.obj$pars_natural <- pp
  mm$.__enclos_env__$private$.obj$state
  mm$.__enclos_env__$private$.obj$update(365*50, 0.1)
  mm$.__enclos_env__$private$.obj$state
#  stopifnot(all.equal(0, sum(mm$.__enclos_env__$private$.obj$state[-(1:2)])))
  n1 <- sum(mm$.__enclos_env__$private$.obj$state[c(3:14)])
  p1 <- sum(mm$.__enclos_env__$private$.obj$state[c("Af","Cf")]) / sum(mm$.__enclos_env__$private$.obj$state[c(3:14)]) * 100
  mm$.__enclos_env__$private$.obj$update(365*1, 0.1)
  n2 <- sum(mm$.__enclos_env__$private$.obj$state[c(3:14)])
  p2 <- sum(mm$.__enclos_env__$private$.obj$state[c("Af","Cf")]) / sum(mm$.__enclos_env__$private$.obj$state[c(3:14)]) * 100
  abs(n1-300) + abs(n2-300) + abs(p1-15) + abs(p2-15)
}

#optim(c(3.0, 0.38), ff, control=list(trace=10))

pars <- c(3.0, 0.38)

#set_pars_natural(parameters);
#double const prev = 0.0205;
#double const N = 257.5;
#m_S.set_sum(N * (1.0-prev));
#m_I.set_sum(N * prev * 0.6);
#m_Af.set_sum(N * prev * 0.3);
#m_Cf.set_sum(N * prev * 0.1);
#m_Z = -N;


mm <- KoalasV2$new()
pp <- mm$.__enclos_env__$private$.obj$parameters
pp["beta"] <- pars[1]
pp["birthrate"] <- pars[2]
pp["passive_intervention_rate"] <- 0.5
mm$.__enclos_env__$private$.obj$parameters <- pp

lapply(seq_len(365*5), \(x){
  mm$.__enclos_env__$private$.obj$update(1, 1/24)
  mm$.__enclos_env__$private$.obj$state
}) |>
  bind_rows() ->
  output

output
output |> tail()

output |>
  select(S:Z) |>
  apply(1,sum) |>
  sapply(\(x) all.equal(x, 0)) |>
  stopifnot()

output |>
  mutate(Total = -Z, Fertile = S+V+I+N+R, Infectious=I+Af+Cf+If, Immune=V+R+Vf+Rf, Diseased=Af+Cf) |>
  select(Year, Day, Total:Diseased) |>
  pivot_longer(Total:Diseased) |>
  ggplot(aes(x=Year + Day/365, y=value, col=name)) +
  geom_line() +
  labs(caption=str_c("beta: ", pars[1], ",  birthrate: ", pars[2]), x="Time (years)", y="Number") +
  scale_x_continuous(breaks=0:10) +
  geom_vline(xintercept=1, lty="dashed")


output |>
  select(Year,Day,S,V,I,R,Af,Cf) |>
  mutate(Total = S+V+I+R+Af+Cf, Prev = 1-(S/Total)) |>
  filter(Year==1 | Prev>=0.15) |>
  head()

theme_set(theme_light())
pdf("outputs_5y.pdf", width=14)
output |>
  select(Year,Day,S,I,R,Af,Cf) |>
  mutate(Total = S+I+R+Af+Cf) |>
  pivot_longer(S:Total) |>
  ggplot(aes(x=Year + Day/365, y=value, col=name)) +
  geom_line() +
  labs(caption=str_c("beta: ", pars[1], ",  birthrate: ", pars[2]), x="Time (years)", y="Number") +
  scale_x_continuous(breaks=0:10) +
  geom_vline(xintercept=1, lty="dashed")

output |>
  mutate(Total = S+I+R+Af+Cf, All = (I+Af+Cf)/Total*100, I = I/Total*100, Af = Af/Total*100, Cf = Cf/Total*100) |>
  select(Year,Day,I,Af,Cf,All) |>
  pivot_longer(I:All) |>
  ggplot(aes(x=Year + Day/365, y=value, col=name)) +
  geom_line() +
  labs(caption=str_c("beta: ", pars[1], ",  birthrate: ", pars[2]), x="Time (years)", y="Prevalence") +
  scale_y_continuous(breaks=seq(0,100,by=10), limits = c(0,100)) +
  scale_x_continuous(breaks=0:10) +
  geom_vline(xintercept=1, lty="dashed")
dev.off()

output |>
  select(Year,Day,S,I,R,Af,Cf) |>
  mutate(Total = S+I+R+Af+Cf) |>
  tail()

output |>
  select(Year,Day,S,I,R,Af,Cf) |>
  mutate(Total = S+I+R+Af+Cf) |>
  mutate(across(c(I,R,Af,Cf), \(x) x/(Total-S))) |>
  tail()

output |>
  mutate(Total = S+I+R+Af+Cf, All = (I+Af+Cf)/Total*100, I = I/Total*100, Af = Af/Total*100, Cf = Cf/Total*100) |>
  select(Year,Day,I,Af,Cf,All) |>
  tail()
