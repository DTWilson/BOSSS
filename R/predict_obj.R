# Returns predicted objectives based on the fitted models or the
# deterministic objective function.

predict_obj <- function(design, problem, solution){
  # Return a matrix with n_samp rows, each giving objective values
  # at design, sampling from stochastic objective models.
  # Note - design can be a matrix
  obj_vals <- matrix(rep(NA, nrow(design)*nrow(problem$objectives)), nrow = nrow(design))
  for(i in 1:nrow(problem$objectives)){

    out_i <- problem$objectives[i, "out_i"]
    hyp_i <- problem$objectives[i, "hyp_i"]

    model_index <- which(solution$to_model$out_i == out_i & solution$to_model$hyp_i == hyp_i)

    if(problem$objectives$stoch[i]){
      ## For stochastic objectives, we use the model to predict
      p <- DiceKriging::predict.km(solution$models[[model_index]],
                                   newdata=design[,1:problem$dimen, drop=F], type="UK",
                                   light.return = TRUE)
      f <- p$mean
      obj_vals[,i] <- f*problem$objectives$weight[i]
    } else {
      ## For deterministic objectives, use the user-written function
      ## For now, assume independent of hypothesis
      f <- as.numeric(problem$det_obj(design)[problem$objectives[i, "out_i"]])
      obj_vals[,i] <- f*objectives$weight[i]
    }
  }
  return(obj_vals)
}


# Returns a matrix of predicted weighted objective values after the next
# evaluation to be used
# to calculate expected hypervolume improvement. When objectives
# are deterministic the entries will be constant. If all objectives
# are deterministic, collapse to a single prediction.

# To do: use a single set of standard normal samples for each stochastic
# objective, transforming these to the correct predicted distribution.

predict_next_obj <- function(n_samp, design, problem, solution){
  # Return a matrix with n_samp rows, each giving objective values
  # at design, sampling from stochastic objective models.
  # Note - design is a single vector
  obj_vals <- matrix(rep(NA, n_samp*nrow(problem$objectives)), nrow = n_samp)
  for(i in 1:nrow(problem$objectives)){

    out_i <- problem$objectives[i, "out_i"]
    hyp_i <- problem$objectives[i, "hyp_i"]

    model_index <- which(solution$to_model$out_i == out_i & solution$to_model$hyp_i == hyp_i)

    if(problem$objectives$stoch[i]){
      ## For stochastic objectives, we use the reinterpolated model to predict
      p <- DiceKriging::predict.km(solution$models_reint[[model_index]],
                                   newdata=design[,1:problem$dimen, drop=F], type="UK",
                                   light.return = TRUE)
      f <- stats::rnorm(n_samp, p$mean, p$sd)
      obj_vals[,i] <- f*problem$objectives$weight[i]
    } else {
      ## For deterministic objectives, use the user-written function
      ## For now, assume independent of hypothesis
      f <- as.numeric(problem$det_obj(design)[problem$objectives[i, "out_i"]])
      obj_vals[,i] <- f*objectives$weight[i]
    }
  }
  return(obj_vals)
}
