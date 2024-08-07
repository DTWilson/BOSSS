
#' Perform one iteration of Bayesian optimisation
#'
#' @param solution current BOSSS solution.
#' @param problem BOSSS problem.
#' @param N number of simulations to use when computing Monte Carlo estimates.
#' @param design optional vector in the design space to be evaluated. If left
#' NULL, an optimal design will be sought.
#'
#' @return An updated BOSSS solution object.
#'
#' @export
iterate <- function(solution, problem, N, design = NULL) {

  if(is.null(design)) {
    opt <- RcppDE::DEoptim(ehi_infill,
                           lower = problem$design_space$lower,
                           upper = problem$design_space$upper,
                           control=list(trace=FALSE, itermax=100, reltol=1e-1, steptol=50),
                           N = N,
                           problem = problem,
                           solution = solution)

    to_eval <- as.numeric(opt$optim$bestmem)
  } else {
    to_eval <- as.numeric(design)
  }

  solution$DoE <- rbind(solution$DoE, c(to_eval, N))

  y <- MC_estimates(to_eval, hypotheses=problem$hypotheses, N=N, sim=problem$simulation, clust=solution$clust)
  if(!is.null(problem$det_func)) {
    y <- c(y, det_values(to_eval, hypotheses=problem$hypotheses, det_func=problem$det_func))
  }

  n_hyp <- ncol(problem$hypotheses)
  out_dimen <- length(y)/(2*n_hyp)
  for(i in 1:n_hyp){
    for(j in 1:out_dimen){
      s <- i*6 - 6 + j*2 - 1
      solution$results[[i,j]] <- rbind(solution$results[[i,j]], y[s:(s+1)])
    }
  }

  mods <- fit_models(solution$DoE, solution$results, solution$to_model, problem)

  solution$models <- mods[1:nrow(solution$to_model)]
  solution$models_reint <- mods[(nrow(solution$to_model)+1):length(mods)]

  pf_out <- pareto_front(solution, problem)
  solution$p_front <- pf_out[[1]]

  p_set <- cbind(solution$DoE, pf_out[[2]])
  p_set <- p_set[sapply(solution$p_front[,ncol(solution$p_front)], function(y) which(p_set[,ncol(p_set)] == y)), 1:problem$dimen]
  obj_vals <- predict_obj(p_set, problem, solution)
  obj_vals <- t(t(obj_vals)/problem$objectives$weight)
  p_set <- cbind(p_set, obj_vals)
  names(p_set)[(problem$dimen + 1):ncol(p_set)] <- problem$objectives$name
  solution$p_set <- p_set

  solution
}
