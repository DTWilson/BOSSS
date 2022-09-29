predict_obj <- function(design, models, objectives, det_obj, to_model)
{
  obj_vals <- matrix(rep(NA, nrow(design)*nrow(objectives)), nrow = nrow(design))
  for(i in 1:nrow(objectives)){
    if(objectives$stoch[i]){
      model_index <- which(to_model$out_i == objectives[i, "out_i"] &
                             to_model$hyp_i == objectives[i, "hyp_i"])
      print(to_model)
      print(model_index)
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design, type="UK")
      f <- p$mean
      obj_vals[,i] <- f*objectives$weight[i]
    } else {
      ## For deterministic objectives, use the user-written function
      ## For now, assume independent of hypothesis
      f <- as.numeric(det_obj(design)[,objectives[i, "out_i"]])
      obj_vals[,i] <- f*objectives$weight[i]
    }
  }
  return(obj_vals)
}

sample_obj <- function(n_samp, design, models_reint, objectives, det_obj, dim, to_model){
  # Return a matrix with n_samp rows, each giving objective values
  # at design, sampling from stochastic objective models.
  obj_vals <- matrix(rep(NA, n_samp*nrow(objectives)), nrow = n_samp)
  for(i in 1:nrow(objectives)){
    if(objectives$stoch[i]){
      ## For stochastic objectives, we use a reinterpolated model
      model_index <- which(to_model$out_i == objectives[i, "out_i"] &
                             to_model$hyp_i == objectives[i, "hyp_i"])
      p <- DiceKriging::predict.km(models_reint[[model_index]], newdata=design[,1:dim, drop=F], type="UK")
      f <- rnorm(n_samp, p$mean, p$sd)
      obj_vals[,i] <- f*objectives$weight[i]
    } else {
      ## For deterministic objectives, use the user-written function
      ## For now, assume independent of hypothesis
      f <- as.numeric(det_obj(design)[objectives[i, "out_i"]])
      obj_vals[,i] <- f*objectives$weight[i]
    }
  }
  return(obj_vals)
}
