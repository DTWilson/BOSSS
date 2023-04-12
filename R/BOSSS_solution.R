# Constructor
new_BOSSS_solution <- function(DoE, results, models, models_reint, p_front, p_set, to_model){

  sol <- list(DoE = DoE,
              results = results,
              models = models,
              models_reint = models_reint,
              p_front = p_front,
              p_set = p_set,
              to_model = to_model)

  structure(sol,
            class = "BOSSS_solution"
  )
}

#' @export
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
  cat("Initialisation will take approximately", substr(dif, 20, nchar(dif)), "\n")

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
  p_front <- pf_out[[1]]
  cat("Initial solution found\n")

  p_set <- cbind(DoE, pf_out[[2]])
  p_set <- p_set[sapply(p_front[,ncol(p_front)], function(y) which(p_set[,ncol(p_set)] == y)), 1:problem$dimen]
  obj_vals <- predict_obj(p_set, models, problem$objectives, problem$det_obj, to_model)
  obj_vals <- t(t(obj_vals)/problem$objectives$weight)
  p_set <- cbind(p_set, obj_vals)
  names(p_set)[(problem$dimen + 1):ncol(p_set)] <- problem$objectives$name

  sol <- new_BOSSS_solution(DoE, results, models, models_reint, p_front, p_set, to_model)
  sol
}

#' @export
print.BOSSS_solution <- function(x, ...) {
  x$p_set
}

#' @export
plot.BOSSS_solution <- function(x, y, ...) {

  n_obj <- (ncol(x$p_front) - 1)
  df <- data.frame(x$p_set)
  obj_names <- names(df)[(ncol(df) - n_obj + 1):ncol(df)]
  names(df)[(ncol(df) - n_obj + 1):ncol(df)] <- letters[1:n_obj]

  if(n_obj > 3){
    stop(
      "Plotting methods are not currently available for > 3 objectives",
      call. = FALSE
    )
  }
  if(n_obj == 3){
    ggplot2::ggplot(df, ggplot2::aes(x=a, y=b, colour=c)) + ggplot2::geom_point() +
      ggplot2::xlab(obj_names[1]) + ggplot2::ylab(obj_names[2]) +
      #viridis::scale_color_viridis(name=obj_names[3]) +
      ggplot2::scale_colour_gradient(name = obj_names[3],
                            low = "green", high = "blue", na.value = NA) +
      ggplot2::theme_minimal()

  } else if(n_obj == 2) {
    ggplot2::ggplot(df, ggplot2::aes(x=a, y=b)) + ggplot2::geom_point() +
      ggplot2::xlab(obj_names[1]) + ggplot2::ylab(obj_names[2]) +
      ggplot2::theme_minimal()

  } else {
    stop(
      "Plotting methods are not currently available for 1 objective",
      call. = FALSE
    )
  }
}
