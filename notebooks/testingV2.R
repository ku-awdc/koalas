library("tidyverse")
library("koalas")

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

pars <- c(3.0, 0.345)

mm <- KoalasV2$new()
pp <- mm$.__enclos_env__$private$.obj$pars_natural
pp["beta"] <- pars[1]
pp["birthrate"] <- pars[2]
mm$.__enclos_env__$private$.obj$pars_natural <- pp

lapply(seq_len(52*20), \(x){
  mm$.__enclos_env__$private$.obj$update(7, 1/24)
  mm$.__enclos_env__$private$.obj$state
}) |>
  bind_rows() ->
  output

theme_set(theme_light())
pdf("outputs.pdf", width=14)
output |>
  select(Year,Day,S,I,R,Af,Cf) |>
  mutate(Total = S+I+R+Af+Cf) |>
  pivot_longer(S:Total) |>
  ggplot(aes(x=Year + Day/365, y=value, col=name)) +
  geom_line() +
  labs(caption=str_c("beta: ", pars[1], ",  birthrate: ", pars[2]), x="Time (years)", y="Number")

output |>
  mutate(Total = S+I+R+Af+Cf, All = (I+Af+Cf)/Total*100, I = I/Total*100, Af = Af/Total*100, Cf = Cf/Total*100) |>
  select(Year,Day,I,Af,Cf,All) |>
  pivot_longer(I:All) |>
  ggplot(aes(x=Year + Day/365, y=value, col=name)) +
  geom_line() +
  labs(caption=str_c("beta: ", pars[1], ",  birthrate: ", pars[2]), x="Time (years)", y="Prevalence") +
  ylim(0,100)
dev.off()

output |>
  select(Year,Day,S,I,R,Af,Cf) |>
  mutate(Total = S+I+R+Af+Cf) |>
  tail()

output |>
  mutate(Total = S+I+R+Af+Cf, All = (I+Af+Cf)/Total*100, I = I/Total*100, Af = Af/Total*100, Cf = Cf/Total*100) |>
  select(Year,Day,I,Af,Cf,All) |>
  tail()
