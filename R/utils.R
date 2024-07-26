# Set up initial space-filling design
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

# Fit surrogate models
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

# Generate MC estimates
MC_estimates <- function(design, hypotheses, N, sim, clust = NULL)
{
  # Run the simulation under each hypothesis and store the results
  results <- NULL
  hypotheses <- t(hypotheses)
  for(i in 1:nrow(hypotheses)) {

    #n.cores <- parallel::detectCores()
    #clust <- parallel::makeCluster(n.cores)
    #parallel::stopCluster(clust)

    if(is.null(clust)){
      sims <- sapply(1:N,
                     eval.parent(substitute(function(...) sim(design, as.data.frame(hypotheses)[i,]) )),
                     simplify = "array")
    } else {
      sims <- parallel::parSapply(clust, 1:N,
                                  eval.parent(substitute(function(...) sim(design, as.data.frame(hypotheses)[i,]) )),
                                  simplify = "array")
    }

    # Make sure output is in matrix form where row = output
    if(is.null(nrow(sims))){
      output_names <- names(sims[1])
    } else {
      output_names <- rownames(sims)
    }
    sims <- matrix(sims, nrow = length(sims)/N)
    for(j in 1:nrow(sims)) {
      # Results are the mean and variance of each of the simulation outputs
      if(all(sims[j,] %in% 0:1)){
        # If the outcome in binary, use a conjugate Beta for estimating mean
        # and variance to avoid getting estimates of 0 with no uncertainty
        # (could try to address with a nugget in the GP models instead?)
        a <- 0.2 + sum(sims[j,] == 1); b <- 0.2 + sum(sims[j,] == 0)
        m <- a/(a+b); v <- (N*m*(1-m)/(0.2+N)^2)*(1/m + 1/(1-m))^2
        #results <- c(results, a/(a + b), a*b/((a+b)^2*(a+b+1)))
        results <- c(results, log(m/(1-m)), v)
      } else {
        results <- c(results, mean(sims[j,]), stats::var(sims[j,])/N)
      }
      # Use the output variable names to name the result columns
      names(results)[(length(results) -1)] <- paste0(output_names[j], "_m_", rownames(hypotheses)[i])
      names(results)[length(results)] <- paste0(output_names[j], "_v_", rownames(hypotheses)[i])
    }
  }
  results
}

det_values <- function(design, hypotheses, det_func) {
  # Evaluate the objective functions under each hypothesis
  results <- NULL
  hypotheses <- t(hypotheses)
  for(i in 1:nrow(hypotheses)) {
    vals <- det_func(design, as.data.frame(hypotheses)[i,])
    # Handle this depending on if there is 1 or more than one output
    for(j in 1:length(vals)) {
      # Deterministic, so variance is 0
      results <- c(results, vals[j], 0)
      # Use the output variable names to name the result columns
      names(results)[(length(results) -1)] <- paste0(names(vals)[j], "_m_", rownames(hypotheses)[i])
      names(results)[length(results)] <- paste0(names(vals)[j], "_v_", rownames(hypotheses)[i])
    }
  }
  results
}



