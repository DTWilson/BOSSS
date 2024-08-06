# Returns predicted objectives based on the fitted models or the
# deterministic objective function.

predict_obj <- function(design, problem, solution){
  # Note - design can be a matrix
  obj_vals <- matrix(rep(NA, nrow(design)*nrow(problem$objectives)), nrow = nrow(design))
  for(i in 1:nrow(problem$objectives)){

    out <- problem$objectives[i, "out"]
    hyp <- problem$objectives[i, "hyp"]


    if(problem$objectives$stoch[i]){
      ## For stochastic objectives, we use the model to predict
      model_index <- which(solution$to_model$out == out & solution$to_model$hyp == hyp)

      p <- DiceKriging::predict.km(solution$models[[model_index]],
                                   newdata=design[,1:problem$dimen, drop=F], type="UK",
                                   light.return = TRUE)
      f <- p$mean
      if(problem$objectives[i, "binary"]) f <- exp(f)/(exp(f) + 1)
      obj_vals[,i] <- f*problem$objectives$weight[i]
    } else {
      ## For deterministic objectives, use the user-written function
      f <- apply(design, 1, problem$det_func, hypothesis = problem$hypotheses[,hyp])[out,]
      obj_vals[,i] <- f*problem$objectives$weight[i]
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

    out <- problem$objectives[i, "out"]
    hyp <- problem$objectives[i, "hyp"]


    if(problem$objectives$stoch[i]){

      model_index <- which(solution$to_model$out == out & solution$to_model$hyp == hyp)

      ## For stochastic objectives, we use the reinterpolated model to predict
      p <- DiceKriging::predict.km(solution$models_reint[[model_index]],
                                   newdata=design[,1:problem$dimen, drop=F], type="UK",
                                   light.return = TRUE)
      f <- stats::rnorm(n_samp, p$mean, p$sd)
      if(problem$objectives[i, "binary"]) f <- exp(f)/(exp(f) + 1)
      obj_vals[,i] <- f*problem$objectives$weight[i]
    } else {
      ## For deterministic objectives
      f <- apply(design, 1, problem$det_func, hypothesis = problem$hypotheses[,hyp])[out,]
      #f <- as.numeric(problem$det_obj(design)[problem$objectives[i, "out_i"]])
      obj_vals[,i] <- f*problem$objectives$weight[i]
    }
  }
  return(obj_vals)
}
