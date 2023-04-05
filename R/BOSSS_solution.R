# Constructor
new_BOSSS_solution <- function(DoE, models, models_reint, p_front, p_set, to_model){

  sol <- list(DoE = DoE,
              models = models,
              models_reint = models_reint,
              p_front = p_front,
              p_set = p_set,
              to_model = to_model)

  structure(sol,
            class = "BOSSS_solution"
  )
}

BOSSS_solution <- function(size, N, problem){
  stopifnot(class(problem) == "BOSSS_problem")

  DoE <- init_DoE(size, problem$design_space)
  DoE$N <- N

  r <-  t(apply(DoE, 1, calc_rates, hypotheses=problem$hypotheses, N=N, sim=problem$simulation))

  # Put results into a (# hyps) x (# outputs) matrix
  n_hyp <- ncol(problem$hypotheses)
  results <- vector(mode = "list", length = n_hyp*problem$out_dimen)
  for(i in 1:n_hyp){
    for(j in 1:problem$out_dimen){
      s <- i*6 - 6 + j*2 - 1
      e <- j + i*3 - 3
      print(c(i,j,s,e))
      results[[e]]  <- r[, s:(s+1)]
    }
  }
  results <- matrix(results, nrow = n_hyp, byrow = TRUE)

  # For now, assume all outputs in all hypotheses are being modelled
  to_model <- data.frame(out_i = rep(1:problem$out_dimen, each = n_hyp),
                         hyp_i = rep(1:n_hyp, problem$out_dimen))

  to_model <- rbind(problem$constraints[,c("out_i", "hyp_i")],
                    problem$objectives[,c("out_i", "hyp_i")])
  to_model <- unique(to_model)

  mods <- fit_models(DoE, results, to_model, problem)
  models <- mods[1:nrow(to_model)]
  models_reint <- mods[(nrow(to_model)+1):length(mods)]

  pf_out <- pareto_front(models, DoE, results, to_model, problem)
  pf <- pf_out[[1]]

  sol <- new_BOSSS_solution(DoE, models, models_reint, pf, pf, to_model)
  sol
}
