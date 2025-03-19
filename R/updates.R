#' Update problem constraints
#'
#' @param problem BOSSS problem
#' @param number index of the constraint(s) to be updated
#' @param name name of the constraint(s) to be updated
#' @param nom new nominal value(s)
#' @param delta new acceptance threshold(s)
#'
#' @return an updated BOSSS problem object.
#' @export
update_constraint <- function(problem, number = NULL, name = NULL,
                              nom = NULL, delta = NULL)
{
  # Checks
  if(!is.null(name)){
    index <- NULL
    for(i in 1:length(name)){
      index <- c(index, which(problem$constraints$name == name))
    }
    if(!is.null(number)){
      if(number != index) stop("Number and name do not match")
    }
  } else {
    if(is.null(number)){
      stop("Please provide either a number OR a name to identify the constraint")
    } else {
      index <- number
    }
  }

  if(!is.null(nom)){
    if(length(nom) != length(index)) stop("Number of nominal values does not match number of constraints")
    problem$constraints$nom[index] <- nom
  }
  if(!is.null(delta)){
    if(length(delta) != length(index)) stop("Number of nominal values does not match number of constraints")
    problem$constraints$delta[index] <- delta
  }

  return(problem)
}

extend_initial <- function(problem, solution, extra_N = 0, extra_points = 0)
{
  if(extra_N == 0 & extra_points == 0){
    stop("Either more simulations or more initial design points, or both, must be specified")
  }

  old_N <- solution$DoE$N[1]
  new_N <- old_N + extra_N

  if(new_N < old_N) stop("new_N must be bigger than current N")

    if(extra_N == 0){
      # Not adding any more sims to existing points
      new_N <- old_N
    } else {
      extra_N <- new_N - old_N
      r_extras <- NULL
      for(i in 1:solution$size){
        r_extra <- MC_estimates(solution$DoE[i,1:problem$dimen], hypotheses=problem$hypotheses, N=extra_N, sim=problem$simulation)
        r_extras <- rbind(r_extras, r_extra)
      }

      # MC_estimate output loops through hyps and then output names
      n_hyp <- ncol(problem$hypotheses)
      out_dimen <- ncol(r_extras)/(2*n_hyp)

      for(i in 1:(ncol(r_extras)/2)){
        hyp_i <- floor((i-1)/out_dimen) + 1
        hyp <- colnames(problem$hypotheses[, hyp_i, drop = FALSE])

        out_i <- (i - 1) %% out_dimen + 1
        out <- names(problem$simulation())[out_i]

        # Check if output is binary
        bin <- any(problem$constraints$binary[problem$constraints$out == out],
                   problem$objectives$binary[problem$objectives$out == out])

        if(bin){
          # Binary outcome
          old_lo <- solution$results[[hyp, out]][1:solution$size, 1]
          old_a <- (old_N + 0.4)/(1 + exp(-old_lo)) - 0.2
          old_b <- old_N - old_a

          new_lo <- r_extras[,2*i - 1]
          new_a <- (extra_N+0.4)/(1 + exp(-new_lo)) - 0.2
          new_b <- extra_N - new_a

          a <- 0.2 + old_a + new_a; b <- 0.2 + old_b + new_b
          N <- old_N + extra_N
          m <- a/(a+b); v <- (N*m*(1-m)/(0.2+N)^2)*(1/m + 1/(1-m))^2

          solution$results[[hyp, out]][1:solution$size, 1] <- log(m/(1 - m))
          solution$results[[hyp, out]][1:solution$size, 2] <- v
        } else {
          # Continuous outcome
          old_m <- solution$results[[hyp, out]][1:solution$size, 1]
          old_v <- solution$results[[hyp, out]][1:solution$size, 2]

          solution$results[[hyp, out]][1:solution$size, 1] <- (old_N*old_m + extra_N*r_extras[,2*i - 1])/(old_N + extra_N)
          solution$results[[hyp, out]][1:solution$size, 1] <- (old_N^2 * old_v + extra_N^2 * r_extras[,2*i])/(old_N + extra_N)^2
        }
      }
    }

  solution$DoE$N[1:solution$size] <- new_N

  if(extra_points != 0){
    # Add extra_points to the initial Sobol sequence and evaluate using
    # new_N simulations
    extra_DoE <- init_DoE(solution$size + extra_points, problem$design_space)[(solution$size+1):(solution$size + extra_points),]
    extra_DoE$N <- new_N

    # Add to the top of the solutions initial design
    solution$DoE <- rbind(extra_DoE, solution$DoE)

    r_extra <- NULL

    for(i in 1:nrow(extra_DoE)){
      r_extra <- rbind(r_extra, MC_estimates(extra_DoE[i,1:problem$dimen], hypotheses=problem$hypotheses, N=new_N, sim=problem$simulation))
    }

    if(is.null(problem$det_func)) {
      r_extra <- r_sim
    } else {
      r_det <- t(apply(extra_DoE[,1:problem$dimen], 1, det_values, hypotheses=problem$hypotheses, det_func=problem$det_func))
      r_extra <- cbind(r_extra, r_det)
    }

    # Add new results to the (# hyps) x (# outputs) results matrix
    # Put results into a (# hyps) x (# outputs) matrix
    n_hyp <- ncol(problem$hypotheses)
    out_dimen <- ncol(r_extra)/(2*n_hyp)
    for(i in 1:n_hyp){
      for(j in 1:out_dimen){
        s <- i*6 - 6 + j*2 - 1
        e <- j + i*3 - 3
        solution$results[[e]]  <- rbind(r_extra[, s:(s+1)], solution$results[[e]])
      }
    }
  }

  solution <- update_solution(solution, problem)

  return(solution)
}
