predict_obj <- function(design, models, objectives, det_obj, dim, to_model)
{
  obj_vals <- NULL
  i <- 1
  for(i in 1:nrow(objectives)){
    if(objectives$stoch[i]){
      ## For stochastic objectives, we are just taking the lowest quantile and
      ## treating this as deterministic - need to improve
      model_index <- which(to_model$out_i == objectives[i, "out_i"] &
                             to_model$hyp_i == objectives[i, "hyp_i"])
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design[,1:dim, drop=F], type="SK")
      f <- p$mean - qnorm(0.7)*p$sd
      obj_vals <- c(obj_vals, f*objectives$weight)
    } else {
      ## For deterministic objectives, use the user-written function
      ## For now, assume independent of hypothesis
      f <- det_obj(design)[objectives[i, "out_i"]]
      obj_vals <- c(obj_vals, f*objectives$weight[i])
    }
  }
  t(objectives$weight*t(matrix(obj_vals, ncol = nrow(objectives))))
}

sample_obj <- function(n_samp, design, models, objectives, det_obj, to_model){
  # Return a matrix with n_samp rows, each giving objective values
  # at design, sampling from stochastic objective models.
  obj_vals <- matrix(rep(NA, n_samp*nrow(objectives)), nrow = n_samp)
  for(i in 1:nrow(objectives)){
    if(objectives$stoch[i]){
      ## For stochastic objectives, we are not using a reinterpolated model yet
      model_index <- which(to_model$out_i == objectives[i, "out_i"] &
                             to_model$hyp_i == objectives[i, "hyp_i"])
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design[,1:dim, drop=F], type="SK")
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
