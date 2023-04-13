#' Plot one-dimensional model fits
#'
#' @param design Design to centre plots at.
#' @param problem BOSSS problem.
#' @param solution BOSSS solution.
#'
#' @return A list of plots of size (# models) x (# design variables).
#' @export
#'
#'
one_d_plots <- function(design, problem, solution) {

  plots <- vector(mode = "list", length = nrow(solution$to_model)*problem$dimen)
  count <- 1
  for(i in 1:nrow(solution$to_model)) {

    for(j in 1:problem$dimen) {
      # Create the grid to evaluate, varying only dimension j
      to_eval <- matrix(rep(NA, problem$dimen*100), ncol = problem$dimen)
      for(k in 1:problem$dimen) {
        if(k == j) {
          range <- c(problem$design_space[k, "lower"], problem$design_space[k, "upper"])
          to_eval[,k] <- seq(range[1], range[2], length.out = 100)
        } else {
          to_eval[,k] <- rep(as.numeric(design[k]), 100)
        }
      }
      to_eval <- as.data.frame(to_eval)
      names(to_eval) <- problem$design_space$name
      to_eval <- rbind(to_eval, design[1:problem$dimen])

      p <- DiceKriging::predict.km(solution$models[[i]],
                                   newdata=to_eval, type="UK",
                                   light.return = TRUE)
      df <- data.frame(m = p$mean, sd = p$sd, v = to_eval[,j])

      pl <- ggplot2::ggplot(df, ggplot2::aes(v, m)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = m - 1.96*sd, ymax = m + 1.96*sd), alpha = 0.2) +
        ggplot2::geom_line() +
        ggplot2::ylab(paste("Mean outcome", solution$to_model$out_i[i], ", hypothesis", solution$to_model$hyp_i[i])) +
        ggplot2::xlab(names(design)[j]) +
        ggplot2::geom_point(data = df[nrow(df),]) +
        ggplot2::theme_minimal()

      plots[[count]] <- pl
      count <- count + 1
    }
  }
  return(plots)
}

#' Empirical check
#'
#' @param design Design to evaluate.
#' @param problem BOSSS problem.
#' @param solution BOSSS solution.
#' @param N Number of simulations to use when computing Monte Carlo estimates.
#' @param current An optional matrix containing the results of previous
#' check_point calls, to be built upon.
#'
#' @return A matrix with each row corresponds to a model of the BOSSS solution
#'  object, giving the Monte Carlo estimates of the mean and variance
#'  along side the number of simulations used to compute them.
#' @export
#'
#'
check_point <- function(design, problem, solution, N, current = NULL) {

  r <- calc_rates(design, problem$hypotheses, N, problem$simulation)

  if(is.null(current)){
    current <- matrix(rep(0, nrow(solution$to_model)*3), ncol = 3)
  }

  for(i in 1:nrow(solution$to_model)) {
    p <- DiceKriging::predict.km(solution$models[[i]],
                                 newdata=design[1:problem$dimen], type="UK",
                                 light.return = TRUE)
    cat(paste0("Model ", i, " prediction interval: [", round(p$lower95, 3), ", ", round(p$upper95, 3), "]\n"))

    out_i <- solution$to_model$out_i[i]
    hyp_i <- solution$to_model$hyp_i[i]
    index <- hyp_i*6 - 6 + out_i*2 - 1
    emp_mean <- r[index]
    emp_var <- r[index + 1]

    current[i, 1] <- current[i, 1]*(current[i, 3]/(current[i, 3] + N)) + emp_mean*(N/(current[i, 3] + N))
    current[i, 2] <- current[i, 2]*(current[i, 3]/(current[i, 3] + N))^2 + emp_var*(N/(current[i, 3] + N))^2
    current[i, 3] <- current[i, 3] + N

    emp_lower95 <- current[i, 1] - 1.96*sqrt(current[i, 2])
    emp_upper95 <- current[i, 1] + 1.96*sqrt(current[i, 2])

    cat(paste0("Model ", i, " empirical interval: [", round(emp_lower95, 3), ", ", round(emp_upper95, 3), "]\n\n"))
  }
  return(current)
}
