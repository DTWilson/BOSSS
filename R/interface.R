# Some helper functions to create the four data frames which partly specify
# the problem.


#' Create a design space
#'
#' @param name Character vector of design variable names.
#' @param lower Numeric vector of lower limits.
#' @param upper Numeric vector of upper limits.
#'
#' @return A data.frame defining the design space.
#' @export
#'
#'
design_space <- function(name, lower, upper) {
  stopifnot(upper > lower)

  data.frame(name = name,
             lower = lower,
             upper = upper)
}


#' Create a set of hypotheses
#'
#' @param name Character vector of hypothesis names.
#' @param param_matrix Numeric matrix, each column giving the model parameter
#' values under a specific hypothesis.
#'
#' @return A data.frame defining the hypotheses.
#' @export
#'
#'
hypotheses <- function(name, param_matrix) {
  df <- as.data.frame(param_matrix)
  row.names(df) <- name
  names(df) <- paste0("hyp_", 1:ncol(param_matrix))
  df
}


#' Create a set of constraints
#'
#' @param name Character vector of constraint names.
#' @param out_i Numeric vector denoting which simulation output each constraint
#' pertains to.
#' @param hyp_i Numeric vector denoting which hypothesis each constraint
#' pertains to.
#' @param nom Numeric vector of nominal upper limits.
#' @param delta Numeric vector of probabilities.
#' @param stoch Boolean vector denoting if the constraint function is stochastic
#' (TRUE) or deterministic (FALSE).
#'
#' @return A data.frame defining the constraints.
#' @export
#'
#'
constraints <- function(name, out_i, hyp_i, nom, delta, stoch) {

  data.frame(name = name,
             out_i= out_i,
             hyp_i = hyp_i,
             nom = nom,
             delta = delta,
             stoch = stoch)
}


#' Create a set of objectives
#'
#' @param name Character vector of objective names.
#' @param out_i Numeric vector denoting which simulation output each objective
#' pertains to.
#' @param hyp_i Numeric vector denoting which hypothesis each objective
#' pertains to.
#' @param weight Numeric vector of weights assigned to each objective.
#' @param stoch Boolean vector denoting if the objective function is stochastic
#' (TRUE) or deterministic (FALSE).
#'
#' @return A data.frame defining the objectives.
#' @export
#'
#'
objectives <- function(name, out_i, hyp_i, weight, stoch) {

  data.frame(name = name,
             out_i= out_i,
             hyp_i = hyp_i,
             weight = weight,
             stoch = stoch)
}
