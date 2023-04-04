# Constructor
new_BOSSS_problem <- function(sim_trial, design_space, hypotheses,
                              constraints, objectives, det_obj){
  # Check types
  stopifnot(is.function(sim_trial))
  stopifnot(is.data.frame(design_space))
  stopifnot(is.data.frame(hypotheses))
  stopifnot(is.data.frame(constraints))
  stopifnot(is.data.frame(objectives))
  stopifnot(is.function(det_obj) | is.null(det_obj))

  # Normalise objective weights
  objectives$weight <- objectives$weight/sum(objectives$weight)

  structure(
    sim_trial,
    design_space = design_space,
    hypotheses = hypotheses,
    constraints = constraints,
    objectives = objectives,
    det_obj = det_obj,
    out_dimen = length(sim_trial()),
    dimen = nrow(design_space),
    class = "BOSSS_problem"
  )
}

# Validator
validate_BOSSS_problem <- function(prob) {

  # For example, check objective output numbers exist
  objectives <- attr(prob, "objectives")
  print(objectives$out_i)

  if(max(objectives$out_i) > attr(prob, "out_dimen")) {
    stop(
      "Objectives refer to non-existing simulation output exist",
      call. = FALSE
    )
  }

  prob
}
