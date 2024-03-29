# Constructor
new_BOSSS_problem <- function(sim_trial, design_space, hypotheses,
                              constraints, objectives, det_func = NULL){
  # Check types
  stopifnot(is.function(sim_trial))
  stopifnot(is.data.frame(design_space))
  stopifnot(is.data.frame(hypotheses))
  stopifnot(is.data.frame(constraints))
  stopifnot(is.data.frame(objectives))
  stopifnot(is.function(det_func) | is.null(det_func))

  # Normalise objective weights
  objectives$weight <- objectives$weight/sum(objectives$weight)

  prob <- list(simulation = sim_trial,
               design_space = design_space,
               hypotheses = hypotheses,
               constraints = constraints,
               objectives = objectives,
               det_func = det_func,
               #out_dimen = length(sim_trial()),
               dimen = nrow(design_space))

  structure(prob,
    class = "BOSSS_problem"
  )
}

# Validator
validate_BOSSS_problem <- function(prob) {

  # Do the functions have default arguments?

  # Does the simulation function return named outputs?
  sim_out_names <- names(prob$simulation())
  det_out_names <- NULL
  if(!is.null(prob$det_func)){
    # Does the deterministic function return named outputs?
    det_out_names <- names(prob$det_func())
  }
  out_names <- c(sim_out_names, det_out_names)
  if(!is.character(out_names)){
    stop("Functions do not return named outputs")
  } else {
    # Do the named outputs appear in the constraints and objectives?
    not_in <- !sapply(out_names, function(x) (x %in% prob$constraints$out) | (x %in% prob$objectives$out))
    if(any(not_in)){
      stop("One or more function outputs not not appear in any constraint or objective")
    }
  }

  # Is the design space defined properly?
  if(any( (prob$design_space$upper - prob$design_space$lower) <= 0 )){
    stop("One or more design space ranges are <= 0")
  }

  # Are constraints defined properly?
  if(any(prob$constraints$delta < 0 | prob$constraints$delta > 1)){
    stop("Constraint deltas must lie in [0,1]")
  }
  if(!is.logical(prob$constraints$stoch)){
    stop("Constraint stochasticity indicators must be logical")
  }

  # Are objectives defined properly?
  if(!is.logical(prob$objectives$stoch)){
    stop("Objective stochasticity indicators must be logical")
  }

  # Do all outputs referred to appear in the function outputs?
  all_req_names <- c(prob$objectives$out, prob$constraints$out)
  not_in <- !sapply(all_req_names, function(x) (x %in% out_names))
  if(any(not_in)){
    stop(paste("Outputs", all_req_names[not_in], " referred to in objectives
               or constraints do not arise as named outputs from
               the simulation or deterministic function"))
  }

  # Is there a deterministic function if required?
  if(any( c(prob$constraints$stoch, prob$objectives$stoch) == F)) {
    if(is.null(prob$det_func)) {
      stop("Determinsitic function required by objectives and or constraints,
           but none supplied")
    }
  }

  prob
}


#' Create a BOSSS problem
#'
#' @param sim_trial Function which generates a single (possibly multivariate)
#' Monte Carlo outcome of a design under a hypothesis.
#' @param design_space Data frame constructed via design_space().
#' @param hypotheses Data frame constructed via hypotheses().
#' @param constraints Data frame constructed via constraints().
#' @param objectives Data frame constructed via objectives().
#' @param det_func Optional function which generates deterministic outcomes of a
#' design under a hypothesis.
#'
#' @return An object of class BOSSS_problem.
#' @export
#'
#'
BOSSS_problem <- function(sim_trial, design_space, hypotheses,
                              constraints, objectives, det_func = NULL){

  prob <- new_BOSSS_problem(sim_trial, design_space, hypotheses,
                                        constraints, objectives, det_func)
  validate_BOSSS_problem(prob)
  prob
}
