# Calculate the expected hypervolume improvement for a given point.

ehi_infill <- function(design, N, PS, models, design_space, constraints, objectives, det_obj, out_dim, to_model)
{
  ## Use the expected hypervolume improvement as implemented in GPareto
  ## Here, all objectives are models, and deterministic objective functions
  ## implemented as such via fastfun to have predict and update methods

  dim <- nrow(design_space)

  ## Get objective value of design
  fs <- predict_obj(design, models, objectives, det_obj, dim)

  design <- as.data.frame(matrix(design, ncol = dim))
  names(design) <- design_space$name

  ## Get expected penalisation if we were to evaluate at design,
  ## using the models of constraint functions
  design$exp_pen <- 1
  for(i in 1:nrow(constraints)){
    index <- DoE_index(constraints[i, "out_i"], constraints[i, "hyp_i"],
                       dim, out_dim)
    model_index <- which(to_model$out_i == constraints[i, "out_i"] & to_model$hyp_i == constraints[i, "hyp_i"])
    nom <- constraints[i, "nom"]
    if(constraints[i, "stoch"]){
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design[,1:dim, drop=F], type="SK")
      mc_vars <- 0.25/N #p$mean*(1-p$mean)/N
      pred_q_mean <- p$mean + stats::qnorm(constraints[i, "delta"])*sqrt(mc_vars*(p$sd^2)/(mc_vars+(p$sd^2)))
      pred_q_var <- ((p$sd^2)^2)/(mc_vars+(p$sd^2))
      design$exp_pen <- design$exp_pen*stats::pnorm(nom, pred_q_mean, sqrt(pred_q_var))
    } else {
      pen <- ifelse(DoE[, index] > nom, 0.0000001, 1)
    }
  }

  ## Improvement is quantified by the increase in the dominated hypervolume of
  ## the approximation set if this design was included
  ref <- design_space$up*objectives$weight
  PS2 <- as.matrix(PS[, objectives$name])
  current <- mco::dominatedHypervolume(PS2, ref)
  pos <- apply(fs, 1, function(obj) mco::dominatedHypervolume(as.matrix(rbind(PS2, obj)), ref) )
  #pos <- mco::dominatedHypervolume(as.matrix(rbind(PS2, fs)), ref)
  imp <- (current-pos)*design$exp_pen

  ## Minimising, so keeping negative
  return(imp)
}
