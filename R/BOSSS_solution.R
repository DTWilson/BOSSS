# Constructor
new_BOSSS_solution <- function(DoE, results, models, models_reint, p_front, p_set, to_model, size, problem=problem, clust=NULL){

  sol <- list(DoE = DoE,
              results = results,
              models = models,
              models_reint = models_reint,
              p_front = p_front,
              p_set = p_set,
              to_model = to_model,
              size = size,
              problem = problem,
              clust = clust)

  structure(sol,
            class = "BOSSS_solution"
  )
}

#' Create an initial solution to a BOSSS problem
#'
#' @param size Number of points used in initialisation.
#' @param N Number of simulations used in MC estimation.
#' @param problem BOSSS problem to be solved.
#'
#' @returns An object of class BOSSS_solution.
#' @export
#'
BOSSS_solution <- function(size, N, problem){
  stopifnot(class(problem) == "BOSSS_problem")

  # Initialise the solution by setting up the initial DoE, evaluating those points,
  # fitting models and extracting the Pareto front and set
  DoE <- init_DoE(size, problem$design_space)
  DoE$N <- N

  #n.cores <- parallel::detectCores()
  #clust <- parallel::makeCluster(n.cores)

  # Get a rough estimate of how long initialisation will take
  cat("\nChecking simulation speed...\n")
  t <- Sys.time()
  r_sim <- try({
    MC_estimates(DoE[1,1:problem$dimen], hypotheses=problem$hypotheses, N=N, sim=problem$simulation)
  }, silent = TRUE)

  if (inherits(r_sim, "try-error")) {
    r_sim
    stop(cat("Simulation threw an error at design ", as.numeric(DoE[1,1:problem$dimen])))
  }


  if(nrow(DoE) > 1) {

    #dif <- utils::capture.output((Sys.time() - t)*(size-1)/n.cores)
    dif <- utils::capture.output((Sys.time() - t)*(size-1))
    cat("Initialisation will take approximately", substr(dif, 20, nchar(dif)), "\n")

    #r_rest <- t(parallel::parApply(clust, DoE[2:nrow(DoE),1:problem$dimen], 1, MC_estimates, hypotheses=problem$hypotheses, N=N, sim=problem$simulation))
    # Use a loop so we can easily identify which (if any) design variables casue
    # the simulation to throw an errorsol$re

    for(i in 2:nrow(DoE)){
      r_next <- try({
        MC_estimates(DoE[i,1:problem$dimen], hypotheses=problem$hypotheses, N=N, sim=problem$simulation)
      }, silent = TRUE)

      if (inherits(r_next, "try-error")) {
        r_next
        stop(cat("Simulation threw an error at design ", as.numeric(DoE[i,1:problem$dimen])))
      }

      r_sim <- rbind(r_sim, r_next)
    }

  }

  if(is.null(problem$det_func)) {
    r <- r_sim
  } else {

    r_det <- t(apply(DoE[,1:problem$dimen], 1, det_values, hypotheses=problem$hypotheses, det_func=problem$det_func))

    r <- cbind(r_sim, r_det)
  }

  # Put results into a (# hyps) x (# outputs) matrix
  n_hyp <- ncol(problem$hypotheses)
  out_dimen <- ncol(r)/(2*n_hyp)
  results <- vector(mode = "list", length = n_hyp*out_dimen)
  for(i in 1:n_hyp){
    for(j in 1:out_dimen){
      s <- i*6 - 6 + j*2 - 1
      e <- j + i*3 - 3
      results[[e]]  <- r[, s:(s+1)]
    }
  }
  results <- matrix(results, nrow = n_hyp, byrow = TRUE)
  rownames(results) <- names(problem$hypotheses)

  if(is.null(problem$det_func)){
    colnames(results) <- names(problem$simulation())
  } else {
    colnames(results) <- c(names(problem$simulation()), names(problem$det_func()))
  }

  # Find the hypothesis x output combinations which need to be modelled
  # (i.e. those forming stochastic constraints or objectives)
  to_model <- rbind(problem$constraints[problem$constraints$stoch, c("out", "hyp")],
                    problem$objectives[problem$objectives$stoch, c("out", "hyp")])
  to_model <- unique(to_model)

  sol <- new_BOSSS_solution(DoE, results, models = NULL, models_reint = NULL, p_front = NULL, p_set = NULL, to_model, size, problem, clust=NULL)

  sol <- update_solution(sol, problem)
  cat("Solution found\n")

  return(sol)
}

# Given the current set of simulation data, fit all GP models and find the
# Pareto set.
update_solution <- function(solution, problem)
{
  mods <- fit_models(solution$DoE, solution$results, solution$to_model, problem)
  solution$models  <- mods[1:nrow(solution$to_model)]
  solution$models_reint <- mods[(nrow(solution$to_model)+1):length(mods)]

  pf_out <- pareto_front(solution, problem)
  solution$p_front <- pf_out[[1]]

  p_set <- cbind(solution$DoE, pf_out[[2]])
  p_set <- p_set[sapply(solution$p_front[,ncol(solution$p_front)], function(y) which(p_set[,ncol(p_set)] == y)), 1:problem$dimen]
  obj_vals <- predict_obj(p_set, problem, solution)
  obj_vals <- t(t(obj_vals)/problem$objectives$weight)
  solution$p_set <- cbind(p_set, obj_vals)
  names(solution$p_set)[(problem$dimen + 1):ncol(solution$p_set)] <- problem$objectives$name

  return(solution)
}


#' Print the Pareto set of a BOSSS solution
#'
#' @param x BOSSS solution.
#' @param ... No other arguments for this method.
#'
#' @return A data.frame containing the Pareto set and the corresponding
#' (unweighted) objective values.
#' @export
#'
#'
print.BOSSS_solution <- function(x, ...) {
  # First, print the Pareto set - i.e. the design variable values
  cat("Design variables for the Pareto set: \n\n")
  print(as.matrix(x$p_set[,1:x$problem$dimen]))

  dfs <- PS_empirical_ests(x, x$problem)
  obj_df <- dfs[[1]]

  cat("\nCorresponding objective function values...\n\n")
  print(obj_df)

  con_df <- dfs[[2]]

  cat("\n...and constraint function values:\n\n")
  print(con_df)
}

#' Title
#'
#' @param solution BOSSS solution.
#' @param problem BOSSS problem.
#'
#' @returns a list of two dataframes with the objective and constraint
#'  estimates (empirical, not model based) of designs in the Pareto set
#' @export
#'
PS_empirical_ests <- function(solution, problem) {
  sol_nums <- as.numeric(rownames(solution$p_set))

  # Next, the corresponding objective values
  objectives <- problem$objectives
  obj_df <- NULL
  labs <- NULL
  for(i in 1:nrow(objectives)){
    out <- objectives[i, "out"]
    hyp <- objectives[i, "hyp"]
    lab <- paste0(out, ", ", hyp, " (mean)")
    labs <- c(labs, lab)

    r <- solution$results[[hyp, out]]
    obj_df <- cbind(obj_df, r[sol_nums, 1])
    colnames(obj_df) <- labs
    if(objectives[i, "stoch"]) {
      lab <- paste0(out, ", ", hyp, " (var)")
      labs <- c(labs, lab)
      obj_df <- cbind(obj_df, r[sol_nums, 2])
      colnames(obj_df) <- labs
    }
  }
  rownames(obj_df) <- sol_nums

  # Finally, the constraint values
  constraints <- problem$constraints
  con_df <- NULL
  labs <- NULL
  for(i in 1:nrow(constraints)){
    out <- constraints[i, "out"]
    hyp <- constraints[i, "hyp"]
    lab <- paste0(out, ", ", hyp, " (mean)")
    labs <- c(labs, lab)

    r <- solution$results[[hyp, out]]
    con_df <- cbind(con_df, r[sol_nums, 1])
    colnames(con_df) <- labs
    if(constraints[i, "stoch"]) {
      lab <- paste0(out, ", ", hyp, " (var)")
      labs <- c(labs, lab)
      con_df <- cbind(con_df, r[sol_nums, 2])
      colnames(con_df) <- labs
    }
  }
  rownames(con_df) <- sol_nums

  return(list(obj_df, con_df))
}

#' Plot the Pareto front of a BOSSS solution
#'
#' @param x BOSSS solution.
#' @param y Not used.
#' @param ... No other arguments for this method.
#'
#' @return A ggplot object two- and three-dimension Pareto fronts.
#' @export
#'
#' @importFrom rlang .data
#'
plot.BOSSS_solution <- function(x, y, ...) {

  n_obj <- (ncol(x$p_front) - 1)
  df <- data.frame(x$p_set)
  obj_names <- x$problem$objectives$name
  names(df)[(ncol(df) - n_obj + 1):ncol(df)] <- letters[1:n_obj]
  df$label <- rownames(df)

  if(n_obj > 3){
    stop(
      "Plotting methods are not currently available for > 3 objectives",
      call. = FALSE
    )
  }
  if(n_obj == 3){
    p <- ggplot2::ggplot(df, ggplot2::aes(x=.data$a, y=.data$b)) +
      ggplot2::geom_point(ggplot2::aes(colour=.data$c)) +
      ggplot2::geom_text(ggplot2::aes(label = .data$label), hjust=-0.5, vjust=-0.5) +
      ggplot2::xlab(obj_names[1]) + ggplot2::ylab(obj_names[2]) +
      #viridis::scale_color_viridis(name=obj_names[3]) +
      ggplot2::scale_colour_gradient(name = obj_names[3],
                            low = "green", high = "blue", na.value = NA) +
      ggplot2::theme_minimal()

  } else if(n_obj == 2) {
    p <- ggplot2::ggplot(df, ggplot2::aes(x=.data$a, y=.data$b)) + ggplot2::geom_point() +
      ggplot2::xlab(obj_names[1]) + ggplot2::ylab(obj_names[2]) +
      ggplot2::theme_minimal()

  } else {
    stop(
      "Plotting methods are not currently available for 1 objective",
      call. = FALSE
    )
  }
  p
}
