
#' Run a one-dimensional sensitivity analysis
#'
#' @param design Design to be tested
#' @param name name of the parameter to be varied
#' @param hypothesis name of the hypothesis around which the SA will be
#' conducted
#' @param lower lower bound of the parameter
#' @param upper upper bound of the parameter
#' @param problem BOSSS problem
#' @param num_eval Number of points to evaluate
#' @param N Number of MC samples to use in each evaluation
#'
#' @return A matrix of estimated means and their variances for each simulation
#' output over the sensitivity analysis parameter range.
#'
#' @export
sensitivity <- function(design, name = NULL, hypothesis = NULL, lower, upper,
                        problem, num_eval = 20, N = 100) {

    p_vals <- seq(lower, upper, length.out = num_eval)

    # Run simulations for the given design at all these hypotheses
    results <- NULL
    for(j in 1:num_eval) {

      hyp <- problem$hypotheses[,hypothesis,drop=F]
      hyp[rownames(hyp) == name,] <- p_vals[j]

      sims <- replicate(N, problem$simulation(as.numeric(design), as.matrix(hyp)))
      r <- NULL
      n <- NULL
      # Handle this depending on if there is 1 or more than one output
      if(is.matrix(sims)) {
        for(k in 1:nrow(sims)) {
          # Results are the mean and variance of each of the simulation outputs
          r <- c(r, c(mean(sims[k,]), stats::var(sims[k,])/N))
          n <- c(n, c(paste0(rownames(sims)[k], "_m"), paste0(rownames(sims)[k], "_v")))
        }
      } else {
        r <- c(r, mean(sims), stats::var(sims)/N)
        n <- c(n, c(paste0(names(sims)[1], "_m"), paste0(names(sims)[1], "_v")))
      }
      results <- rbind(results, r)
    }
    colnames(results) <- n

    results <- cbind(p_vals, results)
    colnames(results)[1] <- name

    return(results)
}
