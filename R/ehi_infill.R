# Calculate the expected hypervolume improvement for a given point.

ehi_infill <- function(design, N, pf, mods, design_space, constraints, objectives, det_obj, out_dim, to_model)
{
  ## Use the expected hypervolume improvement as implemented in GPareto
  ## Here, all objectives are models, and deterministic objective functions
  ## implemented as such via fastfun to have predict and update methods

  models <- mods[1]
  models_reint <- mods[2]

  dim <- nrow(design_space)

  design <- as.data.frame(matrix(design, ncol = dim))
  names(design) <- design_space$name

  exp_pen <- exp_penalty(design, N, models, constraints, dim, out_dim, to_model)

  # Get EHI by sampling from predictive dists of stochastic objectives and
  # taking average of the improvement in dominated hypervolumes

  n_samp <- 20
  samp_fs <- sample_obj(n_samp, design, models_reint, objectives, det_obj, dim, to_model)

  # Hack - fix to make general for any objectives
  ref <- c(design_space$up, 1)*objectives$weight
  current <- emoa::dominated_hypervolume(t(pf), ref)
  pos <- apply(samp_fs, 1, function(obj) emoa::dominated_hypervolume(t(as.matrix(rbind(pf, obj))), ref) )
  imp <- (current-mean(pos))*exp_pen

  ## Minimising, so keeping negative
  return(imp)
}

exp_penalty <- function(design, N, models, constraints, dim, out_dim, to_model){
  ## Get expected penalisation if we were to evaluate at design,
  ## using the models of constraint functions
  exp_pen <- 1
  for(i in 1:nrow(constraints)){
    # Get index of constraint output in DoE
    index <- DoE_index(constraints[i, "out_i"], constraints[i, "hyp_i"],
                       dim, out_dim)
    nom <- constraints[i, "nom"]

    if(constraints[i, "stoch"]){
      # If constraint is stochastic, get index of the appropriate model
      model_index <- which(to_model$out_i == constraints[i, "out_i"] & to_model$hyp_i == constraints[i, "hyp_i"])

      # predict constraint at design point
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design[,1:dim, drop=F], type="SK")

      # assuming worst case MC error, get mean and varinace of the predicted quantile
      mc_vars <- 0.25/N #p$mean*(1-p$mean)/N
      pred_q_mean <- p$mean + stats::qnorm(constraints[i, "delta"])*sqrt(mc_vars*(p$sd^2)/(mc_vars+(p$sd^2)))
      pred_q_var <- ((p$sd^2)^2)/(mc_vars+(p$sd^2))

      # Incorporate prob of constraint violation into penalisation
      exp_pen <- exp_pen*stats::pnorm(nom, pred_q_mean, sqrt(pred_q_var))
    } else {
      # If constraint is deterministic, penalise by ~0 if violated
      exp_pen <- ifelse(DoE[, index] > nom, 0.0000001, 1)
    }
  }
  return(exp_pen)
}

