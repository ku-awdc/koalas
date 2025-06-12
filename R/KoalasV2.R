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
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @importFrom checkmate qassert assert_number matchArg
#' @importFrom methods new
#' @importFrom rlang .data
#' @importFrom checkmate qassert assert_number assert_date
#' @importFrom purrr list_simplify accumulate
#' @importFrom forcats fct
#' @importFrom ggh4x scale_y_facet
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
    #' @param start_date the date corresponding to day 0 of the simulation (only used for outputs)
    #'
    #' @return A new within-group model object
    initialize = function(num = 3L, num_V = num, num_I = num, num_N = num, num_R = num, num_A = num, parameters = list(), state = list(), start_date = "2022-01-01"){

      qassert(num_V, "X1(0,)")
      qassert(num_I, "X1(0,)")
      qassert(num_N, "X1(0,)")
      qassert(num_R, "X1(0,)")
      qassert(num_A, "X1(0,)")

      private$.alpha <- c(V=num_V,I=num_I,N=num_N,R=num_R,A=num_A)

      if(is.character(start_date) || is.POSIXct(start_date)) start_date <- as.Date(start_date)
      assert_date(start_date, any.missing=FALSE, len=1L)
      private$.start_date <- start_date

      ## Note: initialise with default parameters/state so they can be overridden later
      if(all(private$.alpha==1L)){
        private$.obj <- new(KoalaGroupD1, private$.alpha, private$default_parameters(), private$default_state())
      }else if(all(private$.alpha==3L)){
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
    #' @param lifespan_acute #11:  Disease-related mortality rate (replacement for #10) - default is calibrated so that 25% die before entering Cf
    #' @param lifespan_chronic #11:  Disease-related mortality rate (replacement for #10)
    #' @param relative_fecundity #12:  Relative fecundity of diseased animals - NOTE: we ignore males as only females are important for reproduction
    #' @param sensitivity sensitivity of lab test to detect shedding in I(f), Af, Cf
    #' @param specificity specificity of lab test to not detect shedding in other compartments
    #' @param cure_prob_N proportion of N(f) animals with a completed treatment course that are cured of infection (go to R(f))
    #' @param cure_prob_I proportion of I(f) animals with a completed treatment course that are cured of infection (go to R(f))
    #' @param cure_prob_A proportion of Af animals with a completed treatment course that are cured of infection (go to Rf)
    #' @param cure_prob_C proportion of Cf animals with a completed treatment course that are cured of infection (go to Rf)
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
      lifespan_acute,
      lifespan_chronic,
      relative_fecundity,
      sensitivity,
      specificity,
      cure_prob_N,
      cure_prob_I,
      cure_prob_A,
      cure_prob_C,
      vaccine_efficacy,
      vaccine_booster,
      passive_intervention_rate){

      parameters <- private$.obj$parameters |> as.list()

      dr <- "N1(0,]" # Note: durations of Inf mean a rate of 0!
      rt <- "N1[0,)"
      pb <- "N1[0,1]"
      zr <- "N1[0,0]"

      parnames <- c("vacc_immune_duration" = dr,
        "vacc_redshed_duration" = dr,
        "natural_immune_duration" = dr,
        "beta" = rt,
        "subcinical_duration" = dr,
        "subclinical_recover_proportion" = pb,
        "diseased_recover_proportion" = zr,
        "birthrate" = rt,
        "acute_duration" = dr,
        "lifespan_natural" = dr,
        "lifespan_acute" = dr,
        "lifespan_chronic" = dr,
        "relative_fecundity" = pb,
        "sensitivity" = pb,
        "specificity" = pb,
        "cure_prob_N" = pb,
        "cure_prob_I" = pb,
        "cure_prob_A" = pb,
        "cure_prob_C" = pb,
        "vaccine_efficacy" = pb,
        "vaccine_booster" = pb,
        "passive_intervention_rate" = rt)

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
    #' @param Day current day of the simulation (you probably shouldn't change this)
    #' @param SumTx cumulative total number of animals with a successful treatment course (you probably shouldn't change this)
    #' @param SumVx cumulative total number of animals vaccinated - this excludes the animals also treated (you probably shouldn't change this)
    #' @param SumRx cumulative total number of animals removed due to voluntary culling and failure to cure (you probably shouldn't change this)
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
        "Day",
        "SumTx",
        "SumVx",
        "SumRx",
        "SumMx")
      statenames <- rep("N1[0,)", length(sn))
      names(statenames) <- sn
      statenames["Day"] <- "X1[0,)"

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
    #' @param record should the state over time be recorded?
    #' @return self, invisibly
    update = function(n_days = 1L, d_time=1/24, record=TRUE){

      qassert(n_days, "X1(0,)")
      qassert(d_time, "N1(0,)")
      qassert(record, "B1")
      if(!record && length(private$.allres)>0L){
        stop("Unable to stop recording - it has already started!")
      }

      private$check_state()

      private$.obj$update(n_days, d_time, record) |>
        lapply(as.list) |>
        lapply(as_tibble) |>
        bind_rows() ->
        newres

      private$check_state()
      if(record) private$.allres <- c(private$.allres, list(newres))

      invisible(self)
    },

    #' @description
    #' Implement a (one-time) active sampling/capture/testing of all animals
    #'
    #' @param proportion the proportion of animals to test/treat/vaccinate
    #' @param cull_positive the proportion of test-positive animals that will be culled
    #' @param cull_acute the proportion of acute diseased animals that will be culled
    #' @param cull_chronic the proportion of chronic diseased animals that will be culled
    active_intervention = function(proportion, cull_positive = 0.0, cull_acute = 0.2, cull_chronic = 0.3){

      qassert(proportion, "N1[0,1]")
      qassert(cull_positive, "N1[0,1]")
      qassert(cull_acute, "N1[0,1]")
      qassert(cull_chronic, "N1[0,1]")

      private$.obj$active_intervention(proportion, cull_positive, cull_acute, cull_chronic)
      private$check_state()

      private$.interventions <- bind_rows(
        private$.interventions,
        tibble(
          Date = self$date,
          Type = "Active",
          PropPositive = proportion,
          PropAcute = proportion,
          PropChronic = proportion,
          CullPositive = cull_positive,
          CullAcute = cull_acute,
          CullChronic = cull_chronic
        ))


      invisible(self)
    },

    #' @description
    #' Implement a (one-time) targeted capture/treatment of diseased animals
    #'
    #' @param proportion the proportion of diseased animals to identify and treat (overridden by prop_acute and/or prop_chronic, if set)
    #' @param prop_acute the proportion of acutely diseased animals to identify and treat (overrides proportion)
    #' @param prop_chronic the proportion of chronic diseased animals to identify and treat (overrides proportion)
    #' @param cull_positive this parameter is ignored (it is only provided for compatibility with the run method)
    #' @param cull_acute the proportion of acute diseased animals that will be culled
    #' @param cull_chronic the proportion of chronic diseased animals that will be culled
    targeted_intervention = function(proportion, prop_acute = proportion, prop_chronic = proportion, cull_positive = 0, cull_acute = 0.2, cull_chronic = 0.3){

      qassert(proportion, "N1[0,1]")
      qassert(cull_positive, "N1[0,1]")
      qassert(cull_acute, "N1[0,1]")
      qassert(cull_chronic, "N1[0,1]")

      private$.obj$targeted_intervention(prop_acute, prop_chronic, cull_acute, cull_chronic)
      private$check_state()

      private$.interventions <- bind_rows(
        private$.interventions,
        tibble(
          Date = self$date,
          Type = "Targeted",
          PropPositive = 0,
          PropAcute = prop_acute,
          PropChronic = prop_chronic,
          CullPositive = 0,
          CullAcute = cull_acute,
          CullChronic = cull_chronic
        ))

      invisible(self)
    },

    #' @description
    #' Set and burn-in the model for standard scenarios (with current parameters)
    #' @param d_time the desired time step (delta time)
    #' @return self, invisibly
    burnin = function(d_time=1/24){

      stopifnot(self$day == 0L)

      prev <- 0.05
      N <- 275
      days <- 175

      state <- do.call(self$set_state, args=private$default_state() |> as.list())
      self$set_state(S=N*(1-prev), I=N*prev)

      private$.start_date <- as.Date("2022-07-01") - days

      ## Do the burnin:
      self$update(n_days = days, d_time=d_time, record=FALSE)

      ## Then get us to 1st October 2025:
      newdays <- as.numeric(as.Date("2025-10-01") - self$date, units="days")
      self$update(n_days = newdays, d_time=d_time, record=TRUE)

      invisible(self)
    },

    #' @description
    #' Run the model for 1 or more years with the standard scenarios
    #' @param years the number of years to run for
    #' @param frequency the number of sampling events per year
    #' @param prop_active the proportion of animals to test/treat/vaccinate with active sampling at each intervention
    #' @param prop_targeted the proportion of diseased animals to treat with targeted interventions at each intervention
    #' @param d_time the desired time step (delta time)
    #' @param ... additional arguments (cull proportions) passed to the active_intervention method
    #' @return self, invisibly
    run = function(years, frequency, prop_active=0, prop_targeted=0, d_time=1/24, ...){

      qassert(years, "X1(0,)")
      qassert(frequency, "X1[0,)")
      qassert(prop_active, "N1[0,1]")
      qassert(prop_targeted, "N1[0,1]")
      private$.all_run_dates <- c(private$.all_run_dates, self$date)

      for(y in seq_len(years)){

        str_c(
          as.numeric(strftime(self$date, format="%Y"))+1,
          "-",
          strftime(self$date, format="%m-%d")
        ) |>
          as.Date() ->
          tdt
        diy <- as.numeric(tdt - self$date, units="days")

        if(frequency == 0L){
          self$update(diy, d_time=d_time)
        }else{
          intvl <- floor(diy / frequency)
          for(s in seq_len(frequency-1L)){
            self$active_intervention(proportion=prop_active, ...)
            self$targeted_intervention(proportion=prop_targeted, ...)
            self$update(intvl, d_time=d_time)
          }
          self$active_intervention(proportion=prop_active, ...)
          self$targeted_intervention(proportion=prop_targeted, ...)
          self$update(diy - (intvl*(frequency-1L)), d_time=d_time)
        }
      }

      invisible(self)
    },

    #' @description
    #' Print method giving an overview of the current state and parameter values
    #' @return self, invisibly
    print = function(){

      cat("Koala model with following state:\n")
      print(self$state)

      invisible(self)
    },

    #' @description
    #' autoplot method for a default plot
    #'
    #' @param show_treatments option to show/hide treatments
    #' @param prev_line an optional dotted line for target prevalence
    #' @param number_line an optional dotted line for a target stable (or starting) population size - NULL means to use the starting population size
    #' @param ymax a named numeric vector of maximum values for the y axis
    #' @param alphas a named numeric vector of alpha values for each subplot
    #' @param colours a named numeric vector of colour values for each subplot/compartment
    #' @param dot_col the colour to use for the dotted prevalence and number lines
    #'
    #' @return a ggplot2 object
    autoplot = function(show_treatments = TRUE, prev_line = 5, number_line = NULL, ymax = c(Treatment=NA_real_, Number=NA_real_, Prevalence=100), alphas = c(Treatment=0.25, Number=0.5, Prevalence=0.25), colours = c(Diseased="#F8766D", Infectious="#F9C945", Healthy="#619CFF", Treatment="forestgreen", Prevalence="black"), dot_col="grey50"){

      if(is.null(number_line)) number_line <- self$results_long |> filter(Compartment=="Total") |> slice(1L) |> pull(Koalas)

      self$results_wide |>
        select("Date":"Rf") |>
        mutate(Total = rowSums(across(-c("Date", "Day")))) ->
        tdata

      tdata |>
        filter(row_number()<=2L | .data$Total > 0) |> # Remove rows with NA prevalence, but keep first 2 rows to make sure the plot is created
        mutate(Subplot = "Prevalence of Infection (%)", Compartment = "Infectious", Sum = .data$I+.data$If+.data$Af+.data$Cf, Percent = .data$Sum/.data$Total * 100) |>
        select("Date", "Subplot", "Compartment", "Value"="Percent") ->
        prevdata

      tdata |>
        mutate(Healthy=.data$S+.data$V+.data$N+.data$R+.data$Sf+.data$Vf+.data$Nf+.data$Rf, Infectious=.data$I+.data$If, Diseased=.data$Af+.data$Cf) |>
        select("Date", "Total", "Healthy", "Infectious", "Diseased") |>
        pivot_longer("Healthy":"Diseased", names_to="Compartment", values_to="Koalas") |>
        mutate(Compartment = fct(.data$Compartment, levels=c("Diseased","Infectious","Healthy"))) |>
        arrange(.data$Date, .data$Compartment) |>
        group_by(Date) |>
        mutate(Cumulative = accumulate(Koalas, `+`), Lower = lag(Cumulative, default=0)) |>
        ungroup() |>
        mutate(Subplot = "Number of Koalas") ->
        totaldata

      self$treatments |>
        filter(.data$Type=="Treated") |>
        mutate(Subplot = "Cumulative Treatments", Compartment = "Treated") |>
        select("Date", "Subplot", "Compartment", "Value"="Cumulative") ->
        treatdata

      lvs <- c(prevdata$Subplot[1], totaldata$Subplot[1], treatdata$Subplot[1])
      prevdata$Subplot <- fct(prevdata$Subplot, levels=lvs)
      totaldata$Subplot <- fct(totaldata$Subplot, levels=lvs)
      treatdata$Subplot <- fct(treatdata$Subplot, levels=lvs)

      pt <- ggplot() +
        geom_ribbon(data=prevdata, mapping=aes(x=.data$Date, ymin=0, ymax=.data$Value), alpha=alphas["Prevalence"], fill=colours["Prevalence"]) +
        geom_ribbon(data=totaldata, mapping=aes(x=.data$Date, ymin=.data$Lower, ymax=.data$Cumulative, fill=.data$Compartment), alpha=alphas["Number"]) +
        geom_line(data=prevdata, mapping=aes(x=.data$Date, y=.data$Value), col=colours["Prevalence"]) +
        geom_line(data=totaldata, mapping=aes(x=.data$Date, y=.data$Cumulative, col=.data$Compartment)) +
        geom_hline(data=tibble(Subplot = prevdata$Subplot[1], Value = prev_line), mapping=aes(yintercept=Value), lty="dotted", col=dot_col) +
        geom_hline(data=tibble(Subplot = totaldata$Subplot[1], Value = number_line), mapping=aes(yintercept=Value), lty="dotted", col=dot_col)

      if(show_treatments){
        pt <- pt +
          geom_ribbon(data=treatdata, mapping=aes(x=.data$Date, ymin=0, ymax=.data$Value), alpha=alphas["Treatment"], fill=colours["Treatment"]) +
        geom_line(data=treatdata, mapping=aes(x=.data$Date, y=.data$Value), col=colours["Treatment"])
      }

      pt <- pt +
        facet_wrap(~.data$Subplot, ncol=1, scales="free_y") +
        xlab(NULL) + ylab(NULL) +
        scale_fill_manual(values=colours) +
        scale_color_manual(values=colours) +
        guides(fill = guide_legend(reverse = TRUE, title=element_blank()), color = "none") +
        theme_light() +
        theme(
          strip.background = element_rect(fill="grey95"),
          strip.text=element_text(color="black")
        )

      pt +
        scale_y_facet(
          .data$Subplot == prevdata$Subplot[1],
          limits = c(0, ymax["Prevalence"])
        ) +
        scale_y_facet(
          .data$Subplot == totaldata$Subplot[1],
          limits = c(0, ymax["Number"])
        ) +
        scale_y_facet(
          .data$Subplot == treatdata$Subplot[1],
          limits = c(0, ymax["Treatment"])
        ) ->
        pt

      return(pt)
    }

  ),

  private = mlist(

    ## Private fields:
    .obj = NULL,
    .alpha = rep(NA_integer_, 5L),
    .start_date = as.Date(NA_character_),
    .all_run_dates = as.Date(character(0)),
    .allres = list(),
    .interventions = tibble(
      Date = as.Date(character(0L)),
      SampleProportion = numeric(0L),
      CullPositive = numeric(0L),
      CullInfected = numeric(0L),
      CullAcute = numeric(0L)
    ),

    default_parameters = function(){
      list(
        vacc_immune_duration = c(1.0, 0.3, 1.5),  #1
        vacc_redshed_duration = c(0.5, 0.1, 1.0), #2 - RELATIVE TO #1
        natural_immune_duration = c(1.0, 1.0, 1.0), #3 - RELATIVE TO #1
        beta = rep(2.75,3), #4
        subcinical_duration = c(0.5, 0.1, 1.0), #5
        subclinical_recover_proportion = c(0.35, NA_real_, NA_real_),  #6
        diseased_recover_proportion = c(0.0, 0.0, 0.0),  #7
        birthrate = rep(0.38,3), #8
        acute_duration = c(0.4, NA_real_, NA_real_), #9 - 99% progress to Dead/Cf within 1 year (with alpha=3)
        lifespan_natural = c(6.0, 3.0, 12.0), #10
        lifespan_chronic = c(4.8, NA_real_, NA_real_), #11
        lifespan_acute = c(4.0, NA_real_, NA_real_), #11 - 25% die before they reach C - i.e. relative to #9 (hand-calibrated for alpha=3 and acute_duration = 0.4)
        relative_fecundity = c(0.0, 0.0, 0.1), #12 - ignoring males

        sensitivity = c(0.95, 0.90, 1.0),
        specificity = c(0.999, 0.99, 1.0),
        cure_prob_N = c(0.9, NA_real_, NA_real_),
        cure_prob_I = c(0.9, NA_real_, NA_real_),
        cure_prob_A = c(0.75, NA_real_, NA_real_),
        cure_prob_C = c(0.6, NA_real_, NA_real_),
        vaccine_efficacy = c(0.5, NA_real_, NA_real_),
        vaccine_booster = c(1, 1, 1), # relative to vaccine efficacy
        passive_intervention_rate = c(0.03, 0.01, 0.05) #-log(1-c(0.02, 0.01, 0.04))
      ) |>
        lapply(\(x) x[1L]) ->
        pars

      ## Relative to #1:
      pars[["vacc_redshed_duration"]] <- pars[["vacc_immune_duration"]] * pars[["vacc_redshed_duration"]]
      pars[["natural_immune_duration"]] <- pars[["vacc_immune_duration"]] * pars[["natural_immune_duration"]]

      ## Relative to vaccine efficacy:
      pars[["vaccine_booster"]] <- pars[["vaccine_booster"]] * pars[["vaccine_efficacy"]]
      stopifnot(pars[["vaccine_booster"]] <= 1)

      unlist(pars)
    },

    default_state = function(){

      prev <- 0.0#0.0205
      N <- 300.0#257.5

      c(
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
        SumTx = 0.0,
        SumVx = 0.0,
        SumRx = 0.0,
        SumMx = 0.0
      )
    },

    check_state = function(){
      state <- private$.obj$state
      stopifnot(state[-14] >= 0, !is.na(state), is.finite(state))
      state <- state[-c(1, 15:18)]
      stopifnot(all.equal(sum(state),0))
    },

    deep_clone = function(name, value){

      if(name==".obj"){
        newobj <- value$clone()
      }else{
        return(value)
      }

    }

  ),

  ## Active binding functions:
  active = mlist(

    #' @field date the current date of the model (read-only)
    date = function(){
      private$.start_date + self$day
    },

    #' @field day the current day number of the model (read-only)
    day = function(){
      private$.obj$vitals["Day"] |> as.numeric()
    },

    #' @field N the total number of animals alive (read-only)
    N = function(){
      private$.obj$vitals["Alive"] |> as.numeric()
    },

    #' @field prevalence the current prevalence (read-only)
    prevalence = function(){
      as.numeric(private$.obj$vitals["Prevalence"]) * 100
    },

    #' @field state a list representing the current state of the model
    state = function(value){
      if(missing(value)){
        state <- as.list(private$.obj$state[c(3:14,1:2,16:19)-1])
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
      if(length(private$.allres)==0L) stop("Model has not been updated")
      private$.allres |>
        bind_rows() |>
        mutate(Date = private$.start_date + .data$Day) |>
        select("Date", "Day", everything())
    },

    #' @field results_long a data frame of results from the model in long format, with aggregated/summarised compartments (read-only)
    results_long = function(){
      self$results_wide |>
        select(.data$Date:.data$Rf) |>
        mutate(Total = rowSums(across(-c("Date", "Day"))), Sum=.data$Total) |>
        mutate(Healthy=.data$S+.data$V+.data$N+.data$R, Infectious=.data$I+.data$If+.data$Af+.data$Cf, Diseased=.data$Af+.data$Cf, Infertile=.data$Sf+.data$Vf+.data$Nf+.data$Rf+.data$If+.data$Diseased, Immune=.data$V+.data$Vf+.data$R+.data$Rf+.data$N+.data$Nf) |>
        select("Date", "Total", "Healthy":"Immune", "Sum") |>
        pivot_longer("Total":"Immune", names_to="Compartment", values_to="Koalas") |>
        mutate(Percent = .data$Koalas / .data$Sum * 100) |>
        select(-"Sum") |>
        mutate(Compartment = fct(.data$Compartment, levels=rev(c("Healthy","Immune","Infectious","Diseased","Infertile","Total"))))
    },

    #' @field interventions a data frame of intervention time points from the model (read-only)
    interventions = function(){
      private$.interventions
    },

    #' @field run_dates a date vector of dates when run was called, which might be useful for e.g. adding dashed lines to plots (read-only)
    run_dates = function(){
      private$.all_run_dates
    },

    #' @field treatments a data frame of cumulative treatments and vaccinations from the model (read-only)
    treatments = function(){
      self$results_wide |>
        select("Date", "SumTx", "SumVx") |>
        pivot_longer(-"Date", names_to="Type", values_to="Cumulative") |>
        mutate(Type = case_match(.data$Type,
          "SumTx" ~ "Treated",
          "SumVx" ~ "Vaccinated"
        )) |>
        mutate(Type = fct(.data$Type, levels=c("Treated","Vaccinated")))
    }

  )
)
