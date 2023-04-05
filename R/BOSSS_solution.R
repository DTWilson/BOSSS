# Constructor
new_BOSSS_solution <- function(DoE, models, models_reint, p_front, p_set_ids, to_model){

  sol <- list(DoE = DoE,
              models = models,
              models_reint = models_reint,
              p_front = p_front,
              p_set_ids = p_set_ids,
              to_model = to_model)

  structure(sol,
            class = "BOSSS_solution"
  )
}

BOSSS_solution <- function(size, N, problem){
  stopifnot(class(problem) == "BOSSS_problem")

  # Initialise the solution by setting up the initial DoE, evaluating those points,
  # fitting models and extracting the Pareto front and set
  DoE <- init_DoE(size, problem$design_space)
  DoE$N <- N

  # Get a rough estimate of how long initialisation will take
  cat("Checking simulation speed...\n")
  t <- Sys.time()
  r <- calc_rates(DoE[1,], hypotheses=problem$hypotheses, N=N, sim=problem$simulation)
  dif <- capture.output((Sys.time() - t)*size)
  cat("Initialisation will take approximately", substr(ch, 20, nchar(dif)), "\n")

  r <-  t(apply(DoE, 1, calc_rates, hypotheses=problem$hypotheses, N=N, sim=problem$simulation))

  # Put results into a (# hyps) x (# outputs) matrix
  n_hyp <- ncol(problem$hypotheses)
  results <- vector(mode = "list", length = n_hyp*problem$out_dimen)
  for(i in 1:n_hyp){
    for(j in 1:problem$out_dimen){
      s <- i*6 - 6 + j*2 - 1
      e <- j + i*3 - 3
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
  cat("Models fitted\n")

  pf_out <- pareto_front(models, DoE, results, to_model, problem)
  cat("Initial Pareto front and set found\n")

  #DoE2 <- cbind(DoE, pf_out[[2]])
  #ps <- DoE2[sapply(pf[,ncol(pf)], function(x) which(DoE2[,ncol(DoE2)] == x)), 1:problem$dimen]

  sol <- new_BOSSS_solution(DoE, models, models_reint, pf_out[[1]], pf_out[[2]], to_model)
  sol
}
