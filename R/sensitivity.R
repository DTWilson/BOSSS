
#' Run a one-dimensional sensitivity analysis
#'
#' @param sa_hypothesis Data frame defining the SA variable and its bounds
#' @param design Design to be tested
#' @param problem BOSSS problem
#' @param solution BOSSS solution
#' @param num_eval Number of points to evaluate
#' @param N Number of MC samples to use in each evaluation
#'
#' @return A matrix of estimated means and their variances for each simulation
#' output over the SA variable range.
#' @export
#'

sensitivity <- function(sa_hypothesis, design, problem, solution, num_eval = 20, N = 100) {

    # Create the grid to evaluate, varying only the SA parameter
    to_eval <- matrix(rep(problem$hypotheses[,sa_hypothesis$hyp], each = num_eval), ncol = nrow(problem$hypotheses))
    to_eval[, sa_hypothesis$sa_param] <- seq(sa_hypothesis$lower, sa_hypothesis$upper, length.out = num_eval)

    # Run simulations for the given design at all these hypotheses
    results <- NULL
    for(j in 1:num_eval) {
      sims <- replicate(N, problem$simulation(design, as.data.frame(to_eval)[j,]))
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

    results <- cbind(to_eval[, sa_hypothesis$sa_param], results)
    colnames(results)[1] <- rownames(problem$hypotheses)[sa_hypothesis$sa_param]

    return(results)
}
