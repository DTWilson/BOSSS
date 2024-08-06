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

  # Are all function arguments given defaults?


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
#'
#' @examples
#'
#'
#'
#' @export
BOSSS_problem <- function(sim_trial, design_space, hypotheses,
                              constraints, objectives, det_func = NULL){

  internal_sim_trial <- reformat_sim(sim_trial, design_space)
  if(is.null(det_func)) {
    internal_det_func <- NULL
  } else {
    internal_det_func <- reformat_det(det_func, design_space)
  }

  # Flag if constraints / objectives are stochastic
  test_out <- sim_trial()
  sim_out_names <- names(test_out)
  objectives$stoch <- objectives$out %in% sim_out_names
  constraints$stoch <- constraints$out %in% sim_out_names

  # Flag if constraints / objectives are binary
  if(!is.null(det_func)){
    test_out_det <- det_func()
  }

  if(!("binary" %in% names(objectives))) objectives <- check_binary(objectives, test_out, test_out_det)
  if(!("binary" %in% names(constraints))) constraints <- check_binary(constraints, test_out, test_out_det)

  prob <- new_BOSSS_problem(internal_sim_trial, design_space, hypotheses,
                                        constraints, objectives, internal_det_func)
  validate_BOSSS_problem(prob)
  prob
}

# Create an internal version of the simulation function which takes
# design and hypotheses vectors as arguments

reformat_sim <- function(sim_trial, design_space){

  arg_num <- length(formals(sim_trial))
  for(i in 1:arg_num){
    if(is.symbol(formals(sim_trial)[[i]])){
      stop("Determinsitic function missing default argument(s)")
    }
  }

  arg_names <- methods::formalArgs(sim_trial)
  defaults <- as.numeric(formals(sim_trial))
  dim <- nrow(design_space)

  int_sim <- function(design, hypothesis = defaults[(dim+1):length(names)]){

    design <- as.numeric(design); hypothesis <- as.numeric(hypothesis)

    args <- as.list(c(design, hypothesis))
    names(args) <- arg_names

    do.call("sim_trial", args)
  }

  formals(int_sim)$design <- defaults[1:dim]
  formals(int_sim)$hypothesis <- defaults[(dim+1):length(arg_names)]

  int_sim
}

reformat_det <- function(det_func, design_space){

  arg_num <- length(formals(det_func))
  for(i in 1:arg_num){
    if(is.symbol(formals(det_func)[[i]])){
      stop("Determinsitic function missing default argument(s)")
    }
  }

  arg_names <- methods::formalArgs(det_func)
  defaults <- as.numeric(formals(det_func))
  dim <- nrow(design_space)

  int_det <- function(design, hypothesis = defaults[(dim+1):length(names)]){

    design <- as.numeric(design); hypothesis <- as.numeric(hypothesis)

    args <- as.list(c(design, hypothesis))
    names(args) <- arg_names

    do.call("det_func", args)
  }

  formals(int_det)$design <- defaults[1:dim]
  formals(int_det)$hypothesis <- defaults[(dim+1):length(arg_names)]

  int_det
}

check_binary <- function(df, test_out, test_out_det){
  for(i in 1:nrow(df)){
    x <- df$out[i]
    if(x %in% names(test_out)){
      df$binary[i] <- is.logical(test_out[names(test_out) == x])
    } else {
      df$binary[i] <- is.logical(test_out_det[names(test_out_det) == x])
    }
  }
  df
}
