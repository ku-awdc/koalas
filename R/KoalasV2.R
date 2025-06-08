library("R6")

mlist <- IPDMR:::mlist
`$` <- IPDMR:::`$`

#' R6/C++ representation of the Koala model V2
#'
#' @description
#' This is the "V2" model that has substantial differences to the Grogan et al model (\code{\link{KoalasV1}})
#'
#' @import R6
#' @import stringr
#' @import dplyr
#' @import tibble
#' @importFrom methods new
#' @importFrom rlang .data
#' @importFrom checkmate qassert assert_number
#' @importFrom purrr list_simplify
#'
#' @export
KoalasV2 <- R6::R6Class("KoalasV2",

  public = mlist(

    #' @description
    #' Create a new single-group koala model
    #'
    #' @param num number of sub-components for all states (unless overridden) - currently must be either 1 or 3
    #' @param num_V number of sub-components for V (and Vf)
    #' @param num_I number of sub-components for I (and If)
    #' @param num_N number of sub-components for N (and Nf)
    #' @param num_R number of sub-components for R (and Rf)
    #' @param num_A number of sub-components for Af
    #' @param state an initialisation state list - see the set_state method for the allowed values
    #' @param parameters a list of parameter values - see the set_parameters method for the allowed values
    #'
    #' @return A new within-group model object
    initialize = function(num = 3L, num_V = num, num_I = num, num_N = num, num_R = num, num_A = num, parameters = list(), state = list()){

      qassert(num_V, "X1(0,)")
      qassert(num_I, "X1(0,)")
      qassert(num_N, "X1(0,)")
      qassert(num_R, "X1(0,)")
      qassert(num_A, "X1(0,)")

      private$.alpha <- c(V=num_V,I=num_I,N=num_N,R=num_R,A=num_A)

      ## Note: initialise with default parameters/state so they can be overridden later
      if(all(private$.alpha==3L)){
        private$.obj <- new(KoalaGroupD1, private$.alpha, private$default_parameters(), private$default_state())
      }else if(all(private$.alpha==1L)){
        private$.obj <- new(KoalaGroupD3, private$.alpha, private$default_parameters(), private$default_state())
      }else{
        ## TODO: enable
        stop("Sub-compartment sizes must currently be either all 3 or all 1")
        private$.obj <- new(KoalaGroupD0, private$.alpha, private$default_parameters(), private$default_state())
      }

      stopifnot(is.list(parameters), is.list(state))
      parameters <- do.call(self$set_parameters, args=parameters)
      state <- do.call(self$set_state, args=state)

      private$check_state()
      private$.allres <- list(self$state |> as_tibble())

      return(self)
    },

    #' @description
    #' Change one or more current parameter values
    #'
    #' @param vacc_immune_duration #1: Average duration of vaccine-related immunity from infection for susceptibles
    #' @param vacc_redshed_duration #2: Average duration of treatment/vaccine-related reduction in shedding for infecteds, relative to #1
    #' @param natural_immune_duration #3:  Average duration of natural immunity following resolved infection , relative to #1
    #' @param beta #4:  Infection rate (frequency dependent)
    #' @param subcinical_duration #5:  Average duration of subclinical infection before progression to “acute” disease OR spontaneous recovery
    #' @param subclinical_recover_proportion #6:  Proportion of animals that will spontaneously recover, rather than progressing to acute disease
    #' @param diseased_recover_proportion #7:  Spontaneous recovery rate for diseased animals  – ASSUMED DOES NOT HAPPEN SO MUST BE ZERO
    #' @param birthrate #8:  Birth rate (assumed not density-dependent, for now)
    #' @param acute_duration #9:  Average duration of acute (increased mortality) phase before progressing to chronic (normal mortality) phase
    #' @param lifespan_natural #10:  Average lifespan of uninfected koalas (assumed not density-dependent, for now)
    #' @param lifespan_diseased #11:  Disease-related mortality rate (replacement for #10) - default is calibrated so that 25% die before entering Cf
    #' @param relative_fecundity #12:  Relative fecundity of diseased animals - NOTE: we ignore males as only females are important for reproduction
    #' @param sensitivity sensitivity of lab test to detect shedding in I(f), Af, Cf
    #' @param specificity specificity of lab test to not detect shedding in other compartments
    #' @param treatment_dest_R proportion of treated I(f)/Af/Cf animals that are cured of infection (go to R(f))
    #' @param treatment_dest_N proportion of treated I(f)/Af/Cf animals that are released as non-shedding (go to N(f))
    #' @param treatment_dest_IAC proportion of treated I(f)/Af/Cf animals that are released while still shedding (possibly due to a false-negative test)
    #' @param treatment_dest_remove proportion of treated I(f)/Af/Cf animals that are removed permanently due to failure to cure
    #' @param vaccine_efficacy proportion of S(f) animals that have effective vaccination i.e. go to V(f)
    #' @param vaccine_booster proportion of already-vaccinated/immune animals that re-start their time in that category due to the “booster effect” i.e. V(f), N(f), and R(f)
    #' @param passive_intervention_rate rate at which animals are brought in for test/treatment/vaccination passively (can be interpreted as the expected number of times per year each animal is brought in)
    #'
    #' @return self, invisibly
    set_parameters = function(
      vacc_immune_duration,
      vacc_redshed_duration,
      natural_immune_duration,
      beta,
      subcinical_duration,
      subclinical_recover_proportion,
      diseased_recover_proportion,
      birthrate,
      acute_duration,
      lifespan_natural,
      lifespan_diseased,
      relative_fecundity,
      sensitivity,
      specificity,
      treatment_dest_R,
      treatment_dest_N,
      treatment_dest_IAC,
      treatment_dest_remove,
      vaccine_efficacy,
      vaccine_booster,
      passive_intervention_rate){

      parameters <- private$.obj$parameters |> as.list()

      rd <- "N1(0,)"
      pb <- "N1[0,1]"
      zr <- "N1[0,0]"

      parnames <- c("vacc_immune_duration" = rd,
        "vacc_redshed_duration" = rd,
        "natural_immune_duration" = rd,
        "beta" = rd,
        "subcinical_duration" = rd,
        "subclinical_recover_proportion" = pb,
        "diseased_recover_proportion" = zr,
        "birthrate" = rd,
        "acute_duration" = rd,
        "lifespan_natural" = rd,
        "lifespan_diseased" = rd,
        "relative_fecundity" = pb,
        "sensitivity" = pb,
        "specificity" = pb,
        "treatment_dest_R" = pb,
        "treatment_dest_N" = pb,
        "treatment_dest_IAC" = pb,
        "treatment_dest_remove" = pb,
        "vaccine_efficacy" = pb,
        "vaccine_booster" = pb,
        "passive_intervention_rate" = rd)

      argnames <- names(formals(self$set_parameters))
      stopifnot(
        length(parnames) == length(parameters),
        length(argnames) == length(parameters),
        names(parnames) %in% names(parameters),
        names(argnames) %in% names(parameters)
      )

      # Set and check
      for(pp in names(parnames)){
        if(!do.call("missing", list(pp))){
          parameters[[pp]] <- get(pp)
          qassert(parameters[[pp]], parnames[pp])
        }
      }

      # Additional check
      sumdst <- parameters[["treatment_dest_R"]] + parameters[["treatment_dest_N"]] + parameters[["treatment_dest_IAC"]] + parameters[["treatment_dest_remove"]]
      if(!all.equal(sumdst, 1.0)) stop("The treatment_dest parameters must sum to 1")

      private$.obj$parameters <- list_simplify(parameters)

      invisible(self)
    },

    #' @description
    #' Change one or more current state values
    #'
    #' @param S number of S
    #' @param V number of V
    #' @param I number of I
    #' @param N number of N
    #' @param R number of R
    #' @param Af number of Af
    #' @param Cf number of Cf
    #' @param Sf number of Sf
    #' @param Vf number of Vf
    #' @param If number of If
    #' @param Nf number of Nf
    #' @param Rf number of Rf
    #' @param Year current year (you probably shouldn't change this)
    #' @param Day current day (you probably shouldn't change this)
    #' @param SumTx cumulative total number of animals treated (you probably shouldn't change this)
    #' @param SumVx cumulative total number of animals vaccinated (you probably shouldn't change this)
    #' @param SumRx cumulative total number of animals removed due to failure to cure (you probably shouldn't change this)
    #' @param SumMx cumulative mortality due to disease excluding SumRx (you probably shouldn't change this)
    #'
    #' @return self, invisibly
    set_state = function(
      S,
      V,
      I,
      N,
      R,
      Af,
      Cf,
      Sf,
      Vf,
      If,
      Nf,
      Rf,
      Year,
      Day,
      SumTx,
      SumVx,
      SumRx,
      SumMx
    ){

      state <- private$.obj$state |> as.list()
      state[["Z"]] <- NULL

      sn <- c("S",
        "V",
        "I",
        "N",
        "R",
        "Af",
        "Cf",
        "Sf",
        "Vf",
        "If",
        "Nf",
        "Rf",
        "Year",
        "Day",
        "SumTx",
        "SumVx",
        "SumRx",
        "SumMx")
      statenames <- rep("N1[0,)", length(sn))
      names(statenames) <- sn
      statenames[c("Year","Day")] <- "X1[0,)"

      argnames <- names(formals(self$set_state))
      stopifnot(
        length(statenames) == length(state),
        length(argnames) == length(state),
        names(statenames) %in% names(state),
        names(argnames) %in% names(state)
      )

      # Set and check
      for(pp in names(statenames)){
        if(!do.call("missing", list(pp))){
          state[[pp]] <- get(pp)
          qassert(state[[pp]], statenames[pp])
        }
      }

      private$.obj$state <- list_simplify(state)

      invisible(self)
    },

    #' @description
    #' Update the model for one or more day
    #' @param n_days the number of days to update for
    #' @param d_time the desired time step (delta time)
    #' @return self, invisibly
    update = function(n_days = 1L, d_time=1/24){

      qassert(n_days, "X1(0,)")
      qassert(d_time, "N1(0,)")

      seq_len(n_days) |>
        lapply(\(x){
          private$.obj$update(1.0, d_time)
          private$check_state()
          as_tibble(self$state)
        }) |>
        bind_rows() ->
        newres

      private$.allres <- c(private$.allres, list(newres))

      invisible(self)
    },

    #' @description
    #' Implement a (one-time) active sampling/capture/testing of all animals
    #' @param number the (maximum) number of animals to test/treat/vaccinate (ignored if proportion is supplied)
    #' @param proportion the proportion of animals to test/treat/vaccinate
    active_intervention = function(number, proportion){

      if(missing(proportion)){
        qassert(number, "N1(0,]")
        proportion <- number / self$N
        if(proportion>1){
          warning("Requested more vaccinations than available animals")
          proportion <- 1
        }
      }

      qassert(proportion, "N1(0,1)")

      private$.obj$active_intervention(proportion)
      private$check_state()

      invisible(self)
    },

    #' @description
    #' Print method giving an overview of the current state and parameter values
    #' @return self, invisibly
    print = function(){

      cat("Koala model with following state:\n")
      print(self$state)

      invisible(self)
    }
  ),

  private = mlist(

    ## Private fields:
    .obj = NULL,
    .alpha = rep(NA_integer_, 5L),
    .allres = list(),

    default_parameters = function(){
      list(
        vacc_immune_duration = c(1.0, 0.3, 1.5),  #1
        vacc_redshed_duration = c(0.5, 0.1, 1.0), #2 - RELATIVE TO #1
        natural_immune_duration = c(1.0, 1.0, 1.0), #3 - RELATIVE TO #1
        beta = rep(3.0,3), #4
        subcinical_duration = c(0.5, 0.1, 1.0), #5
        subclinical_recover_proportion = c(0.05, NA_real_, NA_real_),  #6
        diseased_recover_proportion = c(0.0, 0.0, 0.0),  #7
        birthrate = rep(0.38,3), #8
        acute_duration = c(0.4, NA_real_, NA_real_), #9
        lifespan_natural = c(5.0, 3.0, 12.0), #10
        lifespan_diseased = rep(0.25,3), #11 - 25% die before they reach C - i.e. relative to #9
        relative_fecundity = c(0.0, 0.0, 0.1), #12 - ignoring males

        sensitivity = c(0.95, 0.90, 1.0),
        specificity = c(0.999, 0.99, 1.0),
        treatment_dest_R = rep(0.6, 3),
        treatment_dest_N = rep(0.2, 3),
        treatment_dest_IAC = rep(0, 3),
        treatment_dest_remove = rep(0.2, 3),
        vaccine_efficacy = c(0.5, NA_real_, NA_real_),
        vaccine_booster = c(1, 1, 1), # relative to vaccine efficacy
        passive_intervention_rate = -log(1-c(0.02, 0.01, 0.04))
      ) |>
        lapply(\(x) x[1L]) ->
        pars

      ## Relative to #1:
      pars[["vacc_redshed_duration"]] <- pars[["vacc_immune_duration"]] * pars[["vacc_redshed_duration"]]
      pars[["natural_immune_duration"]] <- pars[["vacc_immune_duration"]] * pars[["natural_immune_duration"]]

      ## lifespan_diseased is calibrated so that a certain % die before moving on, so is relative to acute duration:
      qassert(private$.alpha["A"], "X1(0,)")
      deadline <- qgamma(0.99, private$.alpha["A"], rate=private$.alpha["A"]/pars[["acute_duration"]])
      duration <- optimise(\(x){
        abs(pars[["lifespan_diseased"]] - pgamma(deadline, private$.alpha["A"], rate=private$.alpha["A"]/x))
      }, c(0, 100))$minimum
      # curve(pgamma(x, 3, rate=3/duration), from=0, to=5); abline(v=deadline); abline(h=pars[["lifespan_diseased"]])
      pars[["lifespan_diseased"]] <- duration

      ## Normalise treatment dest against sensitivity
      pars[["treatment_dest_IAC"]] <- 1 - pars[["sensitivity"]]
      for(nn in c("treatment_dest_R","treatment_dest_N","treatment_dest_remove")) pars[[nn]] <- pars[[nn]] - (1 - pars[["sensitivity"]])/3
      stopifnot(
        pars[c("treatment_dest_R","treatment_dest_N","treatment_dest_remove","treatment_dest_IAC")] |> unlist() |> sum() |> all.equal(1)
      )

      ## Relative to vaccine efficacy:
      pars[["vaccine_booster"]] <- pars[["vaccine_booster"]] * pars[["vaccine_efficacy"]]
      stopifnot(pars[["vaccine_booster"]] <= 1)

      unlist(pars)
    },

    default_state = function(){

      prev <- 0.0205
      N <- 257.5

      c(
        Year = 0.0,
        Day = 0.0,
        S = N * (1.0-prev),
        V = 0.0,
        I = N * prev * 0.6,
        N = 0.0,
        R = 0.0,
        Af = N * prev * 0.3,
        Cf = N * prev * 0.1,
        Sf = 0.0,
        Vf = 0.0,
        If = 0.0,
        Nf = 0.0,
        Rf = 0.0,
        Z = -N,
        SumTx = 0.0,
        SumVx = 0.0,
        SumRx = 0.0,
        SumMx = 0.0
      )
    },

    check_state = function(){
      state <- private$.obj$state
      stopifnot(state[-15] >= 0, !is.na(state), is.finite(state), state["Day"]<=365)
      state <- state[-c(1:2, 16:19)]
      stopifnot(all.equal(sum(state),0))
    }

  ),

  ## Active binding functions:
  active = mlist(

    #' @field time the current time point of the model (read-only)
    time = function(){
      private$.obj$vitals[1:2]
    },

    #' @field N the total number of animals alive (read-only)
    N = function(){
      private$.obj$vitals[3] |> as.numeric()
    },

    #' @field state a list representing the current state of the model
    state = function(value){
      if(missing(value)){
        state <- as.list(private$.obj$state[c(3:14,1:2,16:19)])
        return(state)
      }
      stopifnot(is.list(value))
      do.call(self$set_state, args=value)
    },

    #' @field parameters a list of parameter values
    parameters = function(value){
      if(missing(value)) return(private$.obj$parameters |> as.list())
      stopifnot(is.list(value))
      do.call(self$set_parameters, args=value)
    },

    #' @field results_wide a data frame of results from the model in wide format (read-only)
    results_wide = function(){
      private$.allres |>
        bind_rows() |>
        select(.data$Year, .data$Day, everything())
    },

    #' @field results_long a data frame of results from the model in long format, with aggregated/summarised compartments (read-only)
    results_long = function(){
      stop("TODO")
      self$results_wide |>
        identity()
    }

  )
)
