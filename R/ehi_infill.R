# Calculate the expected hypervolume improvement for a given point.

ehi_infill <- function(design, N, problem, solution)
{
  ## Use the expected hypervolume improvement as implemented in GPareto

  n_mod <- nrow(solution$to_model)
  dim <- problem$dimen

  design <- as.data.frame(matrix(design, ncol = dim))
  names(design) <- problem$design_space$name
  design$N <- N

  exp_pen <- exp_penalty(design, problem, solution)

  # Get EHI by sampling from predictive dists of stochastic objectives and
  # taking average of the improvement in dominated hypervolumes
  n_samp <- 50
  samp_fs <- predict_next_obj(n_samp, design, problem, solution)

  # choose ref point as worst objective val in each dimension
  p_front <- solution$p_front[, 1:(ncol(solution$p_front) - 1), drop = FALSE]
  ref <- apply(p_front, 2, max)
  current <- emoa::dominated_hypervolume(t(p_front), ref)
  pos <- apply(samp_fs, 1, function(obj) emoa::dominated_hypervolume(t(as.matrix(rbind(p_front, obj))), ref) )
  imp <- (current-mean(pos))*exp_pen

  ## Minimising, so keeping negative
  return(imp)
}

exp_penalty <- function(design, problem, solution){
  ## Get expected penalisation if we were to evaluate at design,
  ## using the models of constraint functions
  exp_pen <- 1
  for(i in 1:nrow(problem$constraints)){

    out <- problem$constraints[i, "out"]
    hyp <- problem$constraints[i, "hyp"]

    nom <- problem$constraints[i, "nom"]

    if(problem$constraints[i, "stoch"]){

      # If constraint is stochastic, predict constraint at design point

      model_index <- which(solution$to_model$out == out & solution$to_model$hyp == hyp)

      p <- DiceKriging::predict.km(solution$models[[model_index]],
                                   newdata=design[,1:problem$dimen, drop=F],
                                   type="SK",
                                   light.return=TRUE)

      # assuming worst case MC error, get mean and variance of the predicted quantile
      mc_vars <- 0.25/design$N
      pred_q_mean <- p$mean + stats::qnorm(problem$constraints[i, "delta"])*sqrt(mc_vars*(p$sd^2)/(mc_vars+(p$sd^2)))
      pred_q_var <- ((p$sd^2)^2)/(mc_vars+(p$sd^2))

      # Incorporate prob of constraint violation into penalisation
      exp_pen <- exp_pen*stats::pnorm(nom, pred_q_mean, sqrt(pred_q_var))
    } else {
      # If constraint is deterministic, penalise by ~0 if violated

      # Results are a hyp x out matrix
      ests <- solution$results[hyp_i, out_i]

      exp_pen <- ifelse(ests[, 1] > nom, 0.0000001, 1)
    }
  }
  return(exp_pen)
}

