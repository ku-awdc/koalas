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
#' @importFrom stats rmultinom quantile
#' @import R6
#'
#' @export
KoalasR6 <- R6::R6Class("KoalasR6",

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

      ## Progression through SIDRS:
      transmission_rt <- private$get_transmission_rate()
      qassert(transmission_rt, "R1")
      disease_rt <- private$get_disease_rate()
      qassert(disease_rt, "R1")
      recovery_rt <- private$get_recovery_rate()
      qassert(recovery_rt, "R1")
      reversion_rt <- private$get_reversion_rate()
      qassert(reversion_rt, "R1")

      ## Movement to and from V:
      vaccination_rt <- private$get_vaccination_rate()
      qassert(vaccination_rt, "R5")
      waning_rt <- private$get_waning_rate()
      qassert(waning_rt, "R1")

      ##

      mortalities <- private$get_mortality_rates()
      qassert(mortalities, "R5")

      new_S <-

      leave_S <- apply_rates(private$.S,
        c(transmission, private$.vacc, mortalities["S"]),
        d_time,
        update_type = private$.update_type)
      leave_I <- apply_rates(private$.I,
        c(private$.omega*private$.numE, private$.repl),
        d_time,
        update_type = private$.update_type)
      leave_I <- apply_rates(private$.I,
        c(private$.gamma, private$.repl + private$.cull),
        d_time,
        update_type = private$.update_type)
      leave_R <- apply_rates(private$.R,
        private$.delta+private$.repl,
        d_time,
        update_type = private$.update_type)

      ## As for exercise 3B:
      private$.S <- private$.S + leave_R[,1] + sum(leave_E[,2]) + leave_I[,2] - sum(leave_S)
      private$.E <- private$.E + c(leave_S[,1], leave_E[-private$.numE,1]) - apply(leave_E,1,sum)
      private$.I <- private$.I + leave_E[private$.numE,1] - sum(leave_I)
      private$.R <- private$.R + leave_I[,1] + leave_S[,2] - sum(leave_R)

      private$.time <- private$.time + d_time

      ## New safety feature:
      private$check_state()

      invisible(self)
    },

    #' @description
    #' Update the state of the group for several time points
    #' @param add_time the additional time to add to the current time of the model
    #' @param d_time the desired time step (delta time)
    #' @return a data frame of the model state at each (new) time point
    run = function(add_time, d_time){
      c(
        if(self$time==0) list(self$state) else list(),
        lapply(seq(self$time+d_time, self$time+add_time, by=d_time), function(x){
          self$update(d_time)$state
        })
      ) |>
        bind_rows()
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

    ## Private fields - each has a trailing underscore to avoid name clashes:
    .S = 99,
    .E = numeric(), # The length of this is set at initialisation
    .I = 1,
    .R = 0,
    .N = numeric(), # Set by reset_N method
    .update_type = character(), # Set at initialisation
    .numE = integer(), # Set at initialisation
    .transmission_type = character(), # Set at initialisation

    .time = 0,
    .beta = 0.05,
    .omega = 0.05,
    .gamma = 0.025,
    .delta = 0.005,
    .vacc = 0.001,
    .repl = 0.0001,
    .cull = 0.002,
    .trans_external = 0,
    .group_name = NA_character_,

    ## To store state:
    .saved_public = list(),
    .saved_private = list(),

    ## Private methods:
    check_state = function(){
      ## Use the new compartment_rule private utility method (see below):
      qassert(private$.S, private$compartment_rule())
      qassert(private$.E, private$compartment_rule(private$.numE))
      qassert(private$.I, private$compartment_rule())
      qassert(private$.R, private$compartment_rule())
      qassert(private$.N, private$compartment_rule())

      calcN <- private$.S + sum(private$.E) + private$.I + private$.R
      if(private$.update_type=="stochastic"){
        stopifnot(calcN == private$.N)
      }else if(private$.update_type=="deterministic"){
        stopifnot(isTRUE(all.equal(calcN, private$.N)))
      }else{
        stop("Internal logic error")
      }
    },

    get_transmission_rate = function(){
      if(self$transmission_type=="frequency"){
        trans_internal <- private$.beta * private$.I / private$.N
      }else if(self$transmission_type=="density"){
        trans_internal <- private$.beta * private$.I
      }else{
        stop("Unrecognised transmission type: ", self$transmission_type)
      }
      trans_internal + self$trans_external
    },

    reset_N = function(){
      private$.N <- private$.S + sum(private$.E) + private$.I + private$.R
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
    #' @field E total number of exposed individuals - changing this value will cause the provided total to be distributed evenly (deterministic) or randomly (stochastic) between sub-compartments
    E = function(value){
      if(missing(value)) return(sum(private$.E))
      qassert(value, private$compartment_rule())
      ## For E we distribute input randomly/equally over sub-compartments:
      if(private$.update_type=="stochastic"){
        private$.E <- rmultinom(1, value, rep(1/private$.numE,private$.numE))[,1]
      }else if(private$.update_type=="deterministic"){
        private$.E <- rep(value/private$.numE, private$.numE)
      }else{
        stop("Internal logic error")
      }
      private$reset_N()
    },
    #' @field I number of infected animals
    I = function(value){
      if(missing(value)) return(private$.I)
      qassert(value, private$compartment_rule())
      private$.I <- value
      private$reset_N()
    },
    #' @field R number of recovered animals
    R = function(value){
      if(missing(value)) return(private$.R)
      qassert(value, private$compartment_rule())
      private$.R <- value
      private$reset_N()
    },

    #' @field beta the transmission rate parameter per unit time (must be positive)
    beta = function(value){
      if(missing(value)) return(private$.beta)
      assert_number(value, lower=0)
      private$.beta <- value
    },
    #' @field omega the latent progression rate parameter per unit time (must be positive)
    omega = function(value){
      if(missing(value)) return(private$.omega)
      assert_number(value, lower=0)
      private$.omega <- value
    },
    #' @field gamma the recovery rate parameter per unit time (must be positive)
    gamma = function(value){
      if(missing(value)) return(private$.gamma)
      assert_number(value, lower=0)
      private$.gamma <- value
    },
    #' @field delta the reversion rate parameter per unit time (must be positive)
    delta = function(value){
      if(missing(value)) return(private$.delta)
      assert_number(value, lower=0)
      private$.delta <- value
    },
    #' @field vacc the vaccination rate parameter per unit time (must be positive)
    vacc = function(value){
      if(missing(value)) return(private$.vacc)
      assert_number(value, lower=0)
      private$.vacc <- value
    },
    #' @field repl the replacement rate parameter per unit time (must be positive)
    repl = function(value){
      if(missing(value)) return(private$.repl)
      assert_number(value, lower=0)
      private$.repl <- value
    },
    #' @field cull the targeted culling rate parameter per unit time (must be positive)
    cull = function(value){
      if(missing(value)) return(private$.cull)
      assert_number(value, lower=0)
      private$.cull <- value
    },
    #' @field time the current time point of the model (read-only)
    time = function(){
      private$.time
    },

    #' @field N the total number of animals in the group (read-only)
    N = function(){
      private$.N
    },

    #' @field state a data frame representing the current state of the model (read-only)
    state = function(){
      tibble(Time = self$time, S = self$S, E = self$E, I = self$I, R = self$R)
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
      private$.transmission_type <- match.arg(value,
        choices=c("frequency","density"))
    }
  )
)
