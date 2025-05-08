library("R6")

mlist <- IPDMR:::mlist
`$` <- IPDMR:::`$`

#' R6 representation of the Koala model
#'
#' @description
#' This is a within-group SVIDR(D) model class that can either be run on its own
#' or embedded in a between-group model.
#' Note: D is not a real state (and can be negative), it is just used to ensure
#' that the books are balanced including birth/death
#'
#' @import R6
#' @importFrom IPDMR apply_rates
#' @importFrom stats rmultinom quantile
#' @import stringr
#' @importFrom checkmate qassert assert_number
#'
#' @export
KoalasV1 <- R6::R6Class("KoalasV1",

  public = mlist(

    #' @description
    #' Create a new within-group model
    #' @param update_type one of deterministic or stochastic
    #' @param group_name an optional name for the group (will be included in the output state, if provided)
    #' @return A new within-group model object
    initialize = function(update_type = c("deterministic","stochastic"), group_name=NA_character_){
      ## Checking and partial matching:
      update_type <- match.arg(update_type)
      private$.update_type <- update_type

      private$reset_N()

      ## Save group name:
      qassert(group_name, "s1")
      private$.group_name <- group_name

      return(self)
    },

    #' @description
    #' Update the state of the group for a single time point
    #' @param d_time the desired time step (delta time)
    #' @return self, invisibly
    update = function(d_time){

      assert_number(d_time, lower=0)
      private$check_state()

      # Do birth and death separately, first:
      if(is.finite(private$.carrying_capacity)){
        mort_rt <- private$.mortality_natural + ((private$.birth_rate - private$.mortality_natural) * private$.N/private$.carrying_capacity)
      }else{
        mort_rt <- private$.mortality_natural
      }

      ## Birth (leaving Z):
      zrt <- private$.birth_rate * (private$.N + ((private$.relative_fecundity-1.0) * private$.D))
      if(private$.update_type=="deterministic"){
        leave_Z <- zrt*d_time
      }else if(private$.update_type=="stochastic"){
        leave_Z <- rpois(1, zrt*d_time)
      }else{
        stop("Internal logic error")
      }

      leave_S <- apply_rates(private$.S, mort_rt, d_time, update_type = private$.update_type) |> as.numeric()
      leave_I <- apply_rates(private$.I, mort_rt, d_time, update_type = private$.update_type) |> as.numeric()
      leave_D <- apply_rates(private$.D, private$.mortality_disease+mort_rt, d_time, update_type = private$.update_type) |> as.numeric()
      leave_R <- apply_rates(private$.R, mort_rt, d_time, update_type = private$.update_type) |> as.numeric()
      leave_V <- apply_rates(private$.V, mort_rt, d_time, update_type = private$.update_type) |> as.numeric()

      total_deaths <- leave_S+leave_I+leave_D+leave_R+leave_V
      private$.Z <- private$.Z + total_deaths - leave_Z
      private$.S <- private$.S + leave_Z - leave_S
      private$.I <- private$.I - leave_I
      private$.D <- private$.D - leave_D
      private$.R <- private$.R - leave_R
      private$.V <- private$.V - leave_V
      private$.N <- private$.N + leave_Z - total_deaths

      private$check_state()


      # Then do infection dynamics:
      mort_rt <- 0
      leave_Z <- 0

      ## Leaving S (infection or death):
      if((private$.I + private$.D) == 0){
        trans_rt <- private$.trans_external
      }else{
        trans_rt <- private$.trans_external + private$.beta * (private$.I + private$.D) / private$.N
      }
      leave_S <- apply_rates(private$.S,
                             c(trans_rt, mort_rt),
                             d_time,
                             update_type = private$.update_type)

      ## Leaving I (recovery/disease, or death):
      leave_I <- apply_rates(private$.I,
                             c(private$.sigma, mort_rt),
                             d_time,
                             update_type = private$.update_type)
      ## The sigma (->D or ->S) must be separated afterwards:
      stopifnot(dim(leave_I)==c(1,2))
      if(private$.update_type=="deterministic"){
        leave_I <- cbind(leave_I, leave_I[,1] * private$.acute_recovery_prob)
      }else if(private$.update_type=="stochastic"){
        leave_I <- cbind(leave_I, rbinom(nrow(leave_I), leave_I[,1], private$.acute_recovery_prob))
      }else{
        stop("Internal logic error")
      }
      leave_I[,1] <- leave_I[,1] - leave_I[,3]

      ## Leaving D (recovery, disease-related mortality, general mortality):
      leave_D <- apply_rates(private$.D,
                             c(private$.recovery, 0, mort_rt),
                             d_time,
                             update_type = private$.update_type)

      ## Leaving R (waning immunity, mortality):
      leave_R <- apply_rates(private$.R,
                             c(private$.waning_natural, mort_rt),
                             d_time,
                             update_type = private$.update_type)

      ## Leaving V (waning immunity, mortality):
      leave_V <- apply_rates(private$.V,
                             c(private$.waning_vaccine, mort_rt),
                             d_time,
                             update_type = private$.update_type)

      ## Bookkeeping - NB births/deaths will all be zero here!
      total_deaths <- sum(leave_S[,2], leave_I[,2], leave_D[,2:3], leave_R[,2], leave_V[,2])
      private$.Z <- private$.Z + total_deaths - leave_Z
      private$.S <- private$.S + leave_Z + leave_R[,1] + leave_I[,3] + leave_V[,1] - sum(leave_S)
      private$.I <- private$.I + leave_S[,1] - sum(leave_I)
      private$.D <- private$.D + leave_I[,1] - sum(leave_D)
      private$.R <- private$.R + leave_D[,1] - sum(leave_R)
      private$.V <- private$.V - sum(leave_V)
      private$.N <- private$.N + leave_Z - total_deaths

      private$.time <- private$.time + d_time
      private$check_state()

      invisible(self)

    },

    #' @description
    #' Update the state of the group for several time points
    #' @param add_time the additional time to add to the current time of the model
    #' @param d_time the desired time step (delta time)
    #' @param thin thinning parameter (currently ignored)
    #' @return a data frame of the model state at each (new) time point
    run = function(add_time, d_time, thin=1){
      c(
        if(self$time==0) list(self$state) else list(),
        lapply(seq(self$time+d_time, self$time+add_time, by=d_time), function(x){
          self$update(d_time)$state
        })
      ) |>
        bind_rows()
    },

    #' @description
    #' Implement a (one-time) vaccination effort at the current time point
    #' @param number the (maximum) number of animals to vaccinate (ignored if proportion is supplied)
    #' @param proportion the proportion of animals to vaccinate
    #' @param efficacy the efficacy of the vaccine i.e. probability of the animal moving to V (from S, I or D - R and V do not move)
    vaccinate = function(number, proportion, efficacy){

      stopifnot(private$.update_type=="deterministic")
      if(missing(proportion)){
        qassert(number, "N1(0,]")
        proportion <- number / private$.N
        if(proportion>1){
          warning("Requested more vaccinations than available animals")
          proportion <- 1
        }
      }

      qassert(proportion, "N1(0,1)")
      qassert(efficacy, "N1(0,1)")

      for(cc in c(".S",".I",".D")){
        newv <- private[[cc]] * proportion * efficacy
        private[[cc]] <- private[[cc]] - newv
        private$.V <- private$.V + newv
      }

      private$check_state()

      invisible(self)
    },

    #' @description
    #' Implement a (one-time) active sampling/capture/testing of all animals
    #' @param number the (maximum) number of animals to test and remove (ignored if proportion is supplied)
    #' @param proportion the proportion of animals to test and remove
    #' @param sensitivity the sensitivity of the test
    #' @param specificity the specificity of the test
    active_test = function(number, proportion, sensitivity, specificity){

      stopifnot(private$.update_type=="deterministic")
      if(missing(proportion)){
        qassert(number, "N1(0,]")
        proportion <- number / private$.N
        if(proportion>1){
          warning("Requested more vaccinations than available animals")
          proportion <- 1
        }
      }

      qassert(proportion, "N1(0,1)")
      qassert(sensitivity, "N1(0,1)")
      qassert(specificity, "N1(0,1)")

      for(cc in c(".S",".V",".R")){
        testpos <- private[[cc]] * proportion * (1-specificity)
        private[[cc]] <- private[[cc]] - testpos
        private$.N <- private$.N - testpos
        private$.Cfp <- private$.Cfp + testpos
      }
      for(cc in c(".I",".D")){
        testpos <- private[[cc]] * proportion * sensitivity
        private[[cc]] <- private[[cc]] - testpos
        private$.N <- private$.N - testpos
        private$.Ctp <- private$.Ctp + testpos
      }

      private$check_state()

      invisible(self)
    },

    #' @description
    #' Implement a (one-time) passive sampling/capture/testing of diseased animals
    #' @param number the (maximum) number of animals to test and remove (ignored if proportion is supplied)
    #' @param proportion the proportion of animals to test and remove
    #' @param sensitivity the sensitivity of the test
    #' @param specificity the specificity of the test
    #' @param prevalence_other the prevalence of other conditions that resemble clinical disease
    passive_test = function(number, proportion, sensitivity, specificity, prevalence_other){

      stopifnot(private$.update_type=="deterministic")
      if(missing(proportion)){
        qassert(number, "N1(0,]")
        proportion <- number / private$.N
        if(proportion>1){
          warning("Requested more vaccinations than available animals")
          proportion <- 1
        }
      }

      qassert(proportion, "N1(0,1)")
      qassert(sensitivity, "N1(0,1)")
      qassert(specificity, "N1(0,1)")
      qassert(prevalence_other, "N1(0,1)")

      for(cc in c(".S",".V",".R")){
        testpos <- private[[cc]] * prevalence_other * proportion * (1-specificity)
        private[[cc]] <- private[[cc]] - testpos
        private$.N <- private$.N - testpos
        private$.Cfp <- private$.Cfp + testpos
      }
      for(cc in c(".I")){
        testpos <- private[[cc]] * prevalence_other * proportion * sensitivity
        private[[cc]] <- private[[cc]] - testpos
        private$.N <- private$.N - testpos
        private$.Ctp <- private$.Ctp + testpos
      }
      for(cc in c(".D")){
        testpos <- private[[cc]] * proportion * sensitivity
        private[[cc]] <- private[[cc]] - testpos
        private$.N <- private$.N - testpos
        private$.Ctp <- private$.Ctp + testpos
      }

      private$check_state()
      invisible(self)

    },

    #' @description
    #' Print method giving an overview of the current state and parameter values
    #' @return self, invisibly
    print = function(){
      if(is.na(private$.group_name)){
        cat("An SEIR model with ")
      }else{
        cat("An SEIR model with identifier/name '", private$.group_name, "' and ", sep="")
      }
      cat("the following properties:\n\t",
        "S/E/I/R (N) = ", self$S, "/", self$E, "/", self$I, "/", self$R, " (", self$N, ")\n\t",
        "beta/omega/gamma/delta = ", self$beta, "/", self$omega, "/", self$gamma, "/", self$delta, "\n\t",
        "vacc/repl/cull = ", self$vacc, "/", self$repl, "/", self$cull, "\n\t",
        "E compartments = ", private$.numE, "\n\t",
        "I compartments = 1\n\t",
        "R compartments = 1\n\t",
        "external transmission = ", private$.trans_external, "\n\t",
        "update type = ", private$.update_type, "\n\t",
        "transmission type = ", private$.transmission_type, "\n\t",
        "current time = ", self$time, "\n\t",
        sep="")
      invisible(self)

    }
  ),

  private = mlist(

    ## Private fields:
    .Z = 0,
    .S = 0,
    .I = 0,
    .D = 0,
    .R = 0,
    .V = 0,
    .Ctp = 0,
    .Cfp = 0,
    .N = 0,
    .T = 0,

    .update_type = character(), # Set at initialisation
    .transmission_type = character(), # Set at initialisation
    .group_name = NA_character_,

    .time = 0,

    .mortality_natural = 0.25,
    .mortality_disease = 1.5,
    .carrying_capacity = Inf,
    .birth_rate = 0.38,
    .relative_fecundity = 0.5,

    .beta = 4.0,
    .sigma = 1.33,
    .acute_recovery_prob = 0.29,
    .recovery = 0.05,
    .waning_natural = 0.05,
    .waning_vaccine = 0.05,

    .trans_external = 0,

    ## Private methods:
    check_state = function(){

      if(private$.update_type=="stochastic"){
        qassert(private$.Z, "X1")
      }else if(private$.update_type=="deterministic"){
        qassert(private$.Z, "N1")
      }else{
        stop("Internal logic error")
      }
      qassert(private$.Z, private$compartment_rule() |> str_replace(fixed("[0,)"), ""))
      qassert(private$.S, private$compartment_rule())
      qassert(private$.I, private$compartment_rule())
      qassert(private$.D, private$compartment_rule())
      qassert(private$.R, private$compartment_rule())
      qassert(private$.V, private$compartment_rule())
      qassert(private$.Ctp, private$compartment_rule())
      qassert(private$.Cfp, private$compartment_rule())
      qassert(private$.N, private$compartment_rule())
      qassert(private$.T, private$compartment_rule())

      calcN <- private$.S + private$.I + private$.D + private$.R + private$.V
      if(private$.update_type=="stochastic"){
        stopifnot(calcN == private$.N)
      }else if(private$.update_type=="deterministic"){
        #if(!isTRUE(all.equal(calcN, private$.N))) browser()
        stopifnot(isTRUE(all.equal(calcN, private$.N)))
      }else{
        stop("Internal logic error")
      }

      calcT <- calcN + private$.Z + private$.Ctp + private$.Cfp
      if(private$.update_type=="stochastic"){
        stopifnot(calcT == private$.T)
      }else if(private$.update_type=="deterministic"){
        stopifnot(isTRUE(all.equal(calcT, private$.T)))
      }else{
        stop("Internal logic error")
      }

    },

    reset_N = function(){
      private$.N <- private$.S + private$.I + private$.D + private$.R + private$.V
      private$.T <- private$.N + private$.Z + private$.Ctp + private$.Cfp
      ## Run the check_state method to ensure everything is OK:
      private$check_state()
    },

    ## A utility function to get the correct rule for checking S/E/I/R/N:
    compartment_rule = function(length=1){
      qassert(length, "X1(0,)")
      if(private$.update_type=="stochastic"){
        tt <- "X"
      }else if(private$.update_type=="deterministic"){
        tt <- "N"
      }else{
        stop("Internal logic error")
      }
      str_c(tt, length, "[0,)")
    }
  ),

  ## Active binding functions:
  active = mlist(

    #' @field S number of susceptible animals
    S = function(value){
      if(missing(value)) return(private$.S)
      qassert(value, private$compartment_rule())
      private$.S <- value
      private$reset_N()
    },
    #' @field I total number of infected animals
    I = function(value){
      if(missing(value)) return(sum(private$.I))
      qassert(value, private$compartment_rule())
      private$.I <- value
      private$reset_N()
    },
    #' @field D number of diseased animals
    D = function(value){
      if(missing(value)) return(private$.D)
      qassert(value, private$compartment_rule())
      private$.D <- value
      private$reset_N()
    },
    #' @field R number of recovered animals
    R = function(value){
      if(missing(value)) return(private$.R)
      qassert(value, private$compartment_rule())
      private$.R <- value
      private$reset_N()
    },
    #' @field V number of vaccinated animals
    V = function() return(private$.V),
    #' @field Ctp cumulative number of animals removed as true positives on active/passive sampling
    Ctp = function() return(private$.Ctp),
    #' @field Cfp cumulative number of animals removed as false positives on active/passive sampling
    Cfp = function() return(private$.Cfp),

    #' @field mortality_natural mortality rate from natural causes
    mortality_natural = function(value){
      if(missing(value)) return(private$.mortality_natural)
      qassert(value, "N1[0,)")
      private$.mortality_natural <- value
    },
    #' @field mortality_disease mortality rate from disease
    mortality_disease = function(value){
      if(missing(value)) return(private$.mortality_disease)
      qassert(value, "N1[0,)")
      private$.mortality_disease <- value
    },
    #' @field carrying_capacity number of animals supported by the population (can be Inf or NA)
    carrying_capacity = function(value){
      if(missing(value)) return(private$.carrying_capacity)
      qassert(value, "n1(0,]")
      private$.carrying_capacity <- value
    },
    #' @field birth_rate birth rate
    birth_rate = function(value){
      if(missing(value)) return(private$.birth_rate)
      qassert(value, "N1[0,)")
      private$.birth_rate <- value
    },
    #' @field relative_fecundity fecundity of diseased animals relative to other states (S/I/R/V)
    relative_fecundity = function(value){
      if(missing(value)) return(private$.relative_fecundity)
      qassert(value, "N1[0,1]")
      private$.relative_fecundity <- value
    },

    #' @field beta the transmission rate parameter per unit time (must be positive)
    beta = function(value){
      if(missing(value)) return(private$.beta)
      assert_number(value, lower=0)
      private$.beta <- value
    },
    #' @field sigma the progression rate from I to D or S
    sigma = function(value){
      if(missing(value)) return(private$.sigma)
      assert_number(value, lower=0)
      private$.sigma <- value
    },
    #' @field acute_recovery_prob the probability that a progression from I will go to S rather than D
    acute_recovery_prob = function(value){
      if(missing(value)) return(private$.acute_recovery_prob)
      qassert(value, "N1[0,1]")
      private$.acute_recovery_prob <- value
    },
    #' @field recovery the recovery rate parameter per unit time (must be positive)
    recovery = function(value){
      if(missing(value)) return(private$.recovery)
      assert_number(value, lower=0)
      private$.recovery <- value
    },
    #' @field waning_natural the immune waning rate from R (following natural infection)
    waning_natural = function(value){
      if(missing(value)) return(private$.waning_natural)
      assert_number(value, lower=0)
      private$.waning_natural <- value
    },
    #' @field waning_vaccine the immune waning rate from V (following vaccination)
    waning_vaccine = function(value){
      if(missing(value)) return(private$.waning_vaccine)
      assert_number(value, lower=0)
      private$.waning_vaccine <- value
    },

    #' @field time the current time point of the model (read-only)
    time = function(){
      private$.time
    },
    #' @field N the total number of animals alive (read-only)
    N = function(){
      private$.N
    },

    #' @field state a data frame representing the current state of the model (read-only)
    state = function(){
      tibble(Time = self$time, S = self$S, I = self$I, D = self$D, R = self$R,
             V = self$V, N = self$N, Ctp = self$Ctp, Cfp = self$Cfp)
    },

    #' @field trans_external the external transmission parameter
    trans_external = function(value){
      if(missing(value)) return(private$.trans_external)
      qassert(value, "N1[0,)")
      private$.trans_external <- value
    },

    #' @field transmission_type either frequency or density
    transmission_type = function(value){
      if(missing(value)) return(private$.transmission_type)
      stopifnot(self$time == 0)
      private$.transmission_type <- match.arg(value,
        choices=c("frequency","density"))
    }
  )
)
