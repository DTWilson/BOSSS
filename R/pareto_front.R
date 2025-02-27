
pareto_front <- function(solution, problem)
{
  ## Return the objective values of current Pareto optimal solutions,
  ## penalising constrain violations and considering only solutions
  ## where some evaluation has actually happened

  ## Get objective values
  # Add the objective values for all evaluated points, using point estimates
  # from the fitted models
  obj_v <- predict_obj(solution$DoE[,1:problem$dimen], problem, solution)

  ## Penalise constraint violations
  ## Note - some repetition with penalty in EHI so try to consolidate
  exp_pen <- 1
  if(!is.null(problem$constraints)){
    constraints <- problem$constraints
    for(i in 1:nrow(problem$constraints)){
      pen <- check_constraint(i, solution, problem)
      exp_pen <- exp_pen*pen
    }
  }
  obj_v <- obj_v/exp_pen + 1/exp_pen - 1

  # add an ID variable for extracting the PS later
  obj_v <- cbind(obj_v, sapply(1:nrow(obj_v), function(x) which(order(obj_v[, ncol(obj_v)]) == x)))

  ## Get the Pareto front
  pf <- t(emoa::nondominated_points(t(obj_v)))
  return(list(unique(pf), obj_v[,ncol(obj_v)]))
}

#' @export
check_constraint <- function(i, solution, problem)
{
  dimen <- problem$dimen

  out <- problem$constraints[i, "out"]
  hyp <- problem$constraints[i, "hyp"]

  nom <- problem$constraints[i, "nom"]
  if(problem$constraints$binary[i]) nom <- log(nom/(1-nom))

  pen_prob <- NA

  if(problem$constraints$stoch[i]) {

    # Models are in order of to_model
    model_index <- which(solution$to_model$out == out & solution$to_model$hyp == hyp)

    p <- DiceKriging::predict.km(solution$models[[model_index]],
                                 newdata = solution$DoE[,1:dimen, drop=F],
                                 type="UK", light.return = TRUE)
    pred <- p$mean
    # Get probability that constraint is satisfied
    pen_prob <- stats::pnorm(nom, p$mean, p$sd)
    # Penalise based on this probability being below delta
    pen <- ifelse(pen_prob < problem$constraints[i, "delta"], 0.0000001, 1)
  } else {
    pred <- solution$results[[hyp, out]][,1]
    pen <- ifelse(pred > nom, 0.0000001, 1)
  }
  return(matrix(c(pen, pen_prob, pred), ncol=3))
}


