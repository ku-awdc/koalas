#' Assess different intervention strategies
#'
#' @param koala_model a valid \code{\link{KoalasV2}} model
#' @param years the number of years to update the model (passed to the run method of the \code{\link{KoalasV2}} model)
#' @param frequency a vector of annual sampling frequencies to test
#' @param prop_active a vector of the proportion of animals to test/treat/vaccinate with active sampling at each intervention
#' @param prop_targeted a vector of the proportion of diseased animals to treat with targeted interventions at each intervention
#' @param cl passed to \code{\link[pbapply]{pblapply}} - DO NOT CHANGE THIS ON WINDOWS as it won't properly without shared memory forking
#' @param ... additional arguments (cull proportions) passed to the active_intervention method of the \code{\link{KoalasV2}} model
#'
#' @returns a data frame of final prevalence and population size values at the end of each simulation for each combination of frequency and proportion
#'
#' @importFrom pbapply pblapply
#'
#' @export
assess_interventions <- function(koala_model, years, frequency = c(1,2,3,4), prop_active = seq(0, 1, by=0.05), prop_targeted = 0, cl=NULL, ...){

  stopifnot(inherits(koala_model, c("KoalasV2")))
  qassert(years, "X1(0,)")
  qassert(frequency, "X+(0,365)")
  qassert(prop_active, "N+[0,1]")
  qassert(prop_targeted, "N+[0,1]")
  # Let pblapply deal with cl

  expand_grid(
    Frequency = frequency,
    PropActive = prop_active,
    PropTargeted = prop_targeted
  ) |>
    rowwise() |>
    group_split() |>
    pblapply(
      function(x, mod, years, ...){
        m2 <- mod$clone(deep=TRUE)
        m2[["run"]](years = years, frequency = x[["Frequency"]], prop_active = x[["PropActive"]], prop_targeted = x[["PropTargeted"]], ...)
        x |>
          mutate(Prevalence = m2[["prevalence"]], Koalas = m2[["N"]])
      },
      cl = cl,
      mod = koala_model,
      years = years,
      ...
    ) |>
    bind_rows() |>
    mutate(Frequency = .data[["Frequency"]] |> as.character() |> fct())
}
