#' Set up initial space-filling design
#'
#' @param size integer
#' @param design_space data.frame
#'
#' @return data.frame
#' @export
#'
#' @examples
#' design_space <- data.frame(name = c("n", "k"),
#' low = c(100, 10),
#' up = c(500, 100),
#' int = c(TRUE, TRUE))
#'
#' init_DoE(10, design_space)
init_DoE <- function(size, design_space)
{
  dim <- nrow(design_space)
  DoE <- data.frame(randtoolbox::sobol(size, dim))
  names(DoE) <- design_space$name
  for(i in 1:dim){
    DoE[,i] <-  DoE[,i]*(design_space$up[i]-design_space$low[i]) + design_space$low[i]
  }
  DoE[, design_space$int == 1] <- round(DoE[, design_space$int == 1])

  DoE
}

#' Fit surrogate models
#'
#' @param DoE data.frame
#' @param results Matrix of estimates
#' @param to_model data.frame
#' @param problem BOSSS problem
#'
#' @return list of objects of class km
#' @export
#'
#'
fit_models <- function(DoE, results, to_model, problem)
{
  ## To do: change to updating models if already initialised

  dimen <- problem$dimen

  models <- list()
  models_reint <- list()
  sink("NUL")
  for(i in 1:nrow(to_model)){
    r <- results[[to_model[i,2], to_model[i,1]]]
    models <- append(models, DiceKriging::km(~1, design=DoE[,1:dimen], response=r[, 1],
                                noise.var=r[, 2]))
    # reinterpolated model
    models_reint[[i]] <- DiceKriging::km(~1, design=DoE[,1:dimen],
                                     response=DiceKriging::predict.km(models[[i]], newdata = DoE[,1:dimen], type="UK")$mean)
  }
  sink()

  return(c(models, models_reint))
}

#' Generate MC estimates
#'
#' @param design vector
#' @param hypotheses vector
#' @param N integer
#' @param sim function
#'
#' @return vector
#' @export
#'
#' @examples
#' sim_trial <- function(design, hypothesis)
#' {
#'   n <- design[1]; k <- design[2]
#'   mu <- hypothesis[1]
#'
#'   m <- n/k
#'   s_c <- sqrt(0.05 + 0.95/m)
#'   x0 <- rnorm(k, 0, s_c); x1 <- rnorm(k, mu, s_c)
#'   c(t.test(x0, x1)$p.value >= 0.05, n, k)
#' }
#' design <- c(250, 35)
#' hypotheses <- t(c(0.3))
#' calc_rates(design, hypotheses, N = 100, sim = sim_trial)
calc_rates <- function(design, hypotheses, N, sim)
{
  # Run the simulation under each hypothesis and store the results
  results <- NULL
  hypotheses <- t(hypotheses)
  for(i in 1:nrow(hypotheses)){
    sims <- replicate(N, sim(design, as.data.frame(hypotheses)[i,]))
    for(j in 1:nrow(sims)){
      # Results are the mean and variance of each of the simulation outputs
      results <- c(results, mean(sims[j,]), stats::var(sims[j,])/N)
      # Use the output variable names to name the result columns
      names(results)[(length(results) -1)] <- paste0(rownames(sims)[j], "_m_", rownames(hypotheses)[i])
      names(results)[length(results)] <- paste0(rownames(sims)[j], "_v_", rownames(hypotheses)[i])
    }
  }
  results
}

# Now redundant?
extract_outputs <- function(sim_trial)
{
  str_func <- format(sim_trial)
  found_return <- FALSE
  i <- 1
  while(found_return == FALSE & i <= length(str_func)){
    line <- str_func[i]
    r_index <- regexpr(pattern = "return", format(sim_trial)[i])[[1]]
    i <- i + 1
  }
  if(r_index != -1){
    line_split <- strsplit(line, "=")[[1]]
    n_out <- length(line_split) - 1
    out_names <- rep("", n_out)
    for(i in 1:n_out){
      line_chars <- strsplit(line_split[i], "")[[1]]
      done <- FALSE
      j <- length(line_chars)
      while(!done){
        e <- line_chars[j]
        if(e == " "){
          j <- j - 1
        } else if(e == "(" | e == ","){
          done <- TRUE
        } else {
          out_names[i] <- paste0(e, out_names[i])
          j <- j - 1
        }
      }
    }
  }
  return(out_names)
}





