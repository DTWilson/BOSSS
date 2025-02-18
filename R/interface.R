# Some helper functions to create the four data frames which partly specify
# the problem.


#' Create a design space
#'
#' @param name optional character vector of design variable names.
#' @param sim optional simulation function.
#' @param lower numeric vector of lower limits.
#' @param upper numeric vector of upper limits.
#'
#' @return A data.frame defining the design space.
#'
#' @examples
#' design_space(lower = c(10, 3),
#'              upper = c(500, 50),
#'              name = c("n", "k"))
#'
#' @export
design_space <- function(name = NULL, sim = NULL, lower, upper) {

  if(is.null(name)){
    dim <- length(lower)
    if(is.null(sim)){
      stop(
        "name or sim must be specified.",
        call. = FALSE
      )
    }
    name <- methods::formalArgs(sim)[1:dim]
  }

  stopifnot(upper > lower)

  data.frame(name = name,
             lower = lower,
             upper = upper)
}


#' Create a set of hypotheses
#'
#' @param par_name optional character vector of model parameter names.
#' @param sim optional simulation function.
#' @param values numeric matrix, each column giving the model parameter
#' values under a specific hypothesis.
#' @param hyp_names character vector of hypothesis names.
#'
#' @return A data.frame defining the hypotheses.
#'
#' @examples
#' hypotheses(values = matrix(c(0.3, 0.05, 0.95), ncol = 1),
#'            hyp_names = c("alt"),
#'            par_name = c("mu", "var_u", "var_e"))
#'
#' @export
hypotheses <- function(par_name = NULL, sim = NULL, values, hyp_names) {

  if(is.null(par_name)){
    dim <- nrow(values)
    if(is.null(sim)){
      stop(
        "parname or sim must be specified.",
        call. = FALSE
      )
    }
    num_args <- length(methods::formalArgs(sim))
    par_name <- methods::formalArgs(sim)[(num_args - dim + 1):num_args]
  }

  df <- data.frame(v = values)
  row.names(df) <- par_name
  names(df) <- hyp_names
  df
}


#' Create a set of constraints
#'
#' @param name character vector of constraint names.
#' @param out character vector denoting which simulation output each constraint
#' pertains to.
#' @param hyp character vector denoting which hypothesis each constraint
#' pertains to.
#' @param nom numeric vector of nominal upper limits.
#' @param delta numeric vector of probabilities.
#'
#' @return A data.frame defining the constraints.
#'
#' @examples
#' constraints(name = c("tII"),
#'             out = c("s"),
#'             hyp = c("alt"),
#'             nom = c(0.1),
#'             delta = c(0.95))
#'
#' @export
constraints <- function(name, out, hyp, nom, delta) {

  data.frame(name = name,
             out = out,
             hyp = hyp,
             nom = nom,
             delta = delta)
}


#' Create a set of objectives
#'
#' @param name character vector of objective names.
#' @param out character vector denoting which simulation output each objective
#' pertains to.
#' @param hyp character vector denoting which hypothesis each objective
#' pertains to.
#' @param weight numeric vector of weights assigned to each objective.
#' @param binary optional boolean vector denoting if the output of the objective
#' function is binary (TRUE) or continuous (FALSE).
#'
#' @return A data.frame defining the objectives.
#'
#' @examples objectives(name = c("min_n", "min_k"),
#'                      out = c("n", "k"),
#'                      hyp = c("alt", "alt"),
#'                      weight = c(10, 1))
#'
#' @export
objectives <- function(name, out, hyp, weight, binary = NULL) {

  df <- data.frame(name = name,
             out = out,
             hyp = hyp,
             weight = weight)

  if(!is.null(binary)) df$binary = binary

  df
}
