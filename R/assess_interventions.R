#' Assess different intervention strategies
#'
#' @param koala_model a valid \code{\link{KoalasV2}} model
#' @param years the number of years to update the model (passed to the run method of the \code{\link{KoalasV2}} model)
#' @param frequency a vector of annual sampling frequencies to test
#' @param proportion a vector of sampling proportions to test
#' @param cl passed to \code{\link[pbapply]{pblapply}} - DO NOT CHANGE THIS ON WINDOWS as it won't properly without shared memory forking
#' @param ... additional arguments (cull proportions) passed to the active_intervention method of the \code{\link{KoalasV2}} model
#'
#' @returns a data frame of final prevalence and population size values at the end of each simulation for each combination of frequency and proportion
#'
#' @importFrom pbapply pblapply
#'
#' @export
assess_interventions <- function(koala_model, years, frequency = c(1,2,3,4), proportion = seq(0, 1, by=0.05), cl=NULL, ...){

  stopifnot(inherits(koala_model, c("KoalasV2")))
  qassert(years, "X1(0,)")
  qassert(frequency, "X+(0,365)")
  qassert(proportion, "N+[0,1]")

  expand_grid(
    Frequency = frequency,
    Proportion = proportion
  ) |>
    rowwise() |>
    group_split() |>
    pblapply(
      function(x, mod, years, ...){
        m2 <- mod$clone(deep=TRUE)
        m2[["run"]](years = years, sampling_frequency = x[["Frequency"]], proportion = x[["Proportion"]], ...)
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
