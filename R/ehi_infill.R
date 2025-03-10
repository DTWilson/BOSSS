# Calculate the expected hypervolume improvement for a given point.

ehi_infill <- function(design, N, problem, solution)
{
  ## Use the expected hypervolume improvement as implemented in GPareto

  n_mod <- nrow(solution$to_model)
  dim <- problem$dimen

  design <- as.data.frame(matrix(design, ncol = dim))
  names(design) <- problem$design_space$name
  #design$N <- N

  exp_pen <- exp_penalty(design, problem, solution, N)

  # Get EHI by sampling from predictive dists of stochastic objectives and
  # taking average of the improvement in dominated hypervolumes
  n_samp <- 50
  samp_fs <- predict_next_obj(n_samp, design, problem, solution)


  p_front <- solution$p_front[, 1:(ncol(solution$p_front) - 1), drop = FALSE]
  # Note - need a ref point bigger than the marginal worst case of the current
  # P front, otherwise solutions outwith that box will not be explored
  ref <- apply(p_front, 2, max)
  current <- emoa::dominated_hypervolume(t(p_front), ref)
  pos <- apply(samp_fs, 1, function(obj) emoa::dominated_hypervolume(t(as.matrix(rbind(p_front, obj))), ref) )
  imp <- (current-mean(pos))*exp_pen

  ## Minimising, so keeping negative
  return(imp)
}

exp_penalty <- function(design, problem, solution, N){
  ## Get expected penalisation if we were to evaluate at design,
  ## using the models of constraint functions
  exp_pen <- 1
  if(!is.null(problem$constraints)){
    for(i in 1:nrow(problem$constraints)){

      out <- problem$constraints[i, "out"]
      hyp <- problem$constraints[i, "hyp"]

      nom <- problem$constraints[i, "nom"]
      if(problem$constraints$binary[i]) nom <- log(nom/(1-nom))

      if(problem$constraints[i, "stoch"]){

        # If constraint is stochastic, predict constraint at design point

        model_index <- which(solution$to_model$out == out & solution$to_model$hyp == hyp)

        p <- DiceKriging::predict.km(solution$models[[model_index]],
                                     newdata=design[,1:problem$dimen, drop=F],
                                     type="SK",
                                     light.return=TRUE)

        # assuming worst case MC error, get mean and variance of the predicted quantile
        mc_vars <- (N*0.25*(1-0.25)/(0.2+N)^2)*(1/0.25 + 1/(1-0.25))^2
        pred_q_mean <- p$mean + stats::qnorm(problem$constraints[i, "delta"])*sqrt(mc_vars*(p$sd^2)/(mc_vars+(p$sd^2)))
        pred_q_var <- ((p$sd^2)^2)/(mc_vars+(p$sd^2))

        # Incorporate prob of constraint violation into penalisation
        exp_pen <- exp_pen*stats::pnorm(nom, pred_q_mean, sqrt(pred_q_var))
      } else {
        # If constraint is deterministic, penalise by ~0 if violated

        # Results are a hyp x out matrix
        ests <- solution$results[[hyp, out]]
        dif <- (ests[, 1] > nom)*(ests[, 1] - nom)^2

        exp_pen <- exp_pen*(1/(1 + dif))
      }
    }
  }
  return(exp_pen)
}

