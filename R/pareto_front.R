
pareto_front <- function(solution, problem)
{
  ## Return the objective values of current Pareto optimal solutions,
  ## penalising constrain violations and considering only solutions
  ## where some evaluation has actually happened
  dimen <- problem$dimen

  ## Get objective values
  # Add the objective values for all evaluated points, using point estimates
  # from the fitted models
  obj_v <- predict_obj(solution$DoE[,1:problem$dimen], problem, solution)

  ## Penalise constraint violations
  ## Note - some repetition with penalty in EHI so try to consolidate
  constraints <- problem$constraints
  exp_pen <- 1
  for(i in 1:nrow(problem$constraints)){

    out <- problem$constraints[i, "out"]
    hyp <- problem$constraints[i, "hyp"]

    if(problem$constraints$stoch[i]) {

      # Models are in order of to_model
      model_index <- which(solution$to_model$out == out & solution$to_model$hyp == hyp)

      nom <- problem$constraints[i, "nom"]
      nom <- log(nom/(1-nom))

      p <- DiceKriging::predict.km(solution$models[[model_index]],
                                   newdata = solution$DoE[,1:dimen, drop=F],
                                   type="UK", light.return = TRUE)
      # Get probability that constraint is violated
      pen <- stats::pnorm(nom, p$mean, p$sd)
      # Penalise based on this probability being below delta
      pen <- ifelse(pen < constraints[i, "delta"], 0.0000001, 1)
      exp_pen <- exp_pen*pen
    } else {
      pen <- ifelse(solution$results[[hyp, out]][,1] > nom, 0.0000001, 1)
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



