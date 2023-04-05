iterate <- function(solution, problem, N) {

  opt <- RcppDE::DEoptim(ehi_infill,
                         lower = problem$design_space$lower,
                         upper = problem$design_space$upper,
                         control=list(trace=FALSE, itermax=100, reltol=1e-1, steptol=50),
                         N = N,
                         pf = solution$p_front[,1:nrow(problem$objectives)],
                         models = solution$models,
                         models_reint = solution$models_reint,
                         design_space = problem$design_space,
                         constraints = problem$constraints,
                         objectives = problem$objectives,
                         det_obj = problem$det_obj,
                         out_dim = problem$out_dim,
                         to_model = solution$to_model)
  to_eval <- as.numeric(opt$optim$bestmem)

  solution$DoE <- rbind(solution$DoE, c(to_eval, N))

  y <- calc_rates(to_eval, hypotheses=problem$hypotheses, N=N, sim=problem$simulation)

  n_hyp <- ncol(problem$hypotheses)
  for(i in 1:n_hyp){
    for(j in 1:problem$out_dimen){
      s <- i*6 - 6 + j*2 - 1
      solution$results[[i,j]] <- rbind(solution$results[[i,j]], y[s:(s+1)])
    }
  }

  mods <- fit_models(solution$DoE, solution$results, solution$to_model, problem)

  solution$models <- mods[1:nrow(solution$to_model)]
  solution$models_reint <- mods[(nrow(solution$to_model)+1):length(mods)]

  pf_out <- pareto_front(solution$models, solution$DoE, solution$results, solution$to_model, problem)
  solution$p_front <- pf_out[[1]]

  p_set <- cbind(solution$DoE, pf_out[[2]])
  p_set <- p_set[sapply(solution$p_front[,ncol(solution$p_front)], function(y) which(p_set[,ncol(p_set)] == y)), 1:problem$dimen]
  obj_vals <- predict_obj(p_set, solution$models, problem$objectives, problem$det_obj, solution$to_model)
  obj_vals <- t(t(obj_vals)/problem$objectives$weight)
  p_set <- cbind(p_set, obj_vals)
  names(p_set)[(problem$dimen + 1):ncol(p_set)] <- problem$objectives$name
  solution$p_set <- p_set

  solution
}
