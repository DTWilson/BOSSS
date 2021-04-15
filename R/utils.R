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
  DoE[, design_space$int] <- round(DoE[, design_space$int])

  DoE
}

#' Fit surrogate models
#'
#' @param DoE data.frame
#' @param to_model data.frame
#' @param design_space data.frame
#'
#' @return list of objects of class km
#' @export
#'
#' @examples
#' design_space <- data.frame(name = c("n", "k"),
#' low = c(100, 10),
#' up = c(500, 100),
#' int = c(TRUE, TRUE))
#'
#' DoE <- data.frame(n = c(300, 400, 200),
#' k = c(55, 32, 78),
#' a = c(0.11, 0.10, 0.15),
#' b = c(0.000988, 0.000909, 0.00128),
#' N = c(100, 100, 100))
#'
#' to_model <- data.frame(out_i = c(1),
#' hyp_i = c(1))
#'
#' fit_models(DoE, to_model, design_space)
fit_models <- function(DoE, to_model, design_space)
{
  dim <- nrow(design_space)
  out_dim <- 3

  models <- list()
  for(i in 1:nrow(to_model)){
    response_index <- (to_model[i, "hyp_i"] - 1)*out_dim*2 + (2*to_model[i, "out_i"] - 1) + nrow(design_space)
    models <- append(models, DiceKriging::km(~1, design=DoE[1:dim], response=DoE[, response_index],
                                noise.var=DoE[, response_index + 1]))
  }

  models
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
#' hypotheses <- c(0.3)
#' calc_rates(design, hypotheses, N = 100, sim = sim_trial)
calc_rates <- function(design, hypotheses, N, sim)
{
  results <- NULL
  for(i in 1:nrow(as.data.frame(hypotheses))){
    sims <- replicate(N, sim(design, as.data.frame(hypotheses)[i,]))
    for(j in 1:nrow(sims)){
      results <- c(results, mean(sims[j,]), stats::var(sims[j,])/N)
    }
  }
  names(results) <- letters[1:length(results)]
  results
}

DoE_index <- function(out_i, hyp_i, dim, out_dim)
{
  # For an output and hypothesis number, get the column number in the DoE
  (hyp_i - 1)*out_dim*2 + (2*out_i - 1) + dim
}

is_nondom <- function(x, b, objectives)
{
  i <- 1
  obj_names <- as.character(objectives$name)
  while(i <= nrow(objectives) & nrow(b)!= 0 ){
    ## subset b to those non-dominated solutions which are less than or equal to
    ## x in an objective
    b <- b[b[, obj_names[i] ] <= x[[ obj_names[i] ]],]
    i <- i + 1
  }
  ## If b is now empty, then des is non-dominated
  if(nrow(b)==1 | all(apply(b[,obj_names, drop=F], 2, function(x) length(unique(x)) == 1) == TRUE) ) {
    nondom <- TRUE
  } else {
    nondom <- FALSE
  }
  return(nondom)
}

#' Get the Pareto set
#'
#' @param design_space data.frame
#' @param models list of objects of class km
#' @param DoE data.frame
#' @param objectives data.frame
#' @param constraints data.frame
#' @param to_model data.frame
#' @param b data.frame
#'
#' @return data.frame
#' @export
best <- function(design_space, models, DoE, objectives, constraints, to_model, b=NULL)
{
  ## Return the set of current Pareto optimal solutions,
  ## penalising constrain violations and considering only solutions
  ## where some evaluation has actually happened
  sols <- DoE

  dim <- nrow(design_space)
  out_dim <- max(c(objectives$out_i, constraints$out_i))

  ## Get objective values
  for(i in 1:nrow(objectives)){
    index <- DoE_index(objectives[i, "out_i"], objectives[i, "hyp_i"], dim, out_dim)
    if(objectives[i, "stoch"]){
      model_index <- which(to_model$out_i == objectives[i, "out_i"] & to_model$hyp_i == objectives[i, "hyp_i"])
      sols <- cbind(sols, DiceKriging::predict.km(models[[model_index]], newdata=sols[, 1:dim, drop = F], type="SK")$mean*objectives[i, "weight"])
    } else {
      sols <- cbind(sols, DoE[, index]*objectives[i, "weight"])
    }
    names(sols)[ncol(sols)] <- as.character(objectives[i, "name"])
  }

  ## Penalise constraint violations
  sols$exp_pen <- 1
  for(i in 1:nrow(constraints)){
    index <- DoE_index(constraints[i, "out_i"], constraints[i, "hyp_i"], dim, out_dim)
    model_index <- which(to_model$out_i == constraints[i, "out_i"] & to_model$hyp_i == constraints[i, "hyp_i"])
    nom <- constraints[i, "nom"]
    if(constraints[i, "stoch"]){
      p <- DiceKriging::predict.km(models[[model_index]], newdata=sols[,1:dim, drop=F], type="SK")
      pen <- stats::pnorm(nom, p$mean, p$sd)
      pen <- ifelse(pen < constraints[i, "delta"], 0.0000001, 1)
      sols$exp_pen <- sols$exp_pen*pen
    } else {
      pen <- ifelse(DoE[, index] > nom, 0.0000001, 1)
    }
  }
  sols[, objectives[, "name"]] <- sols[, objectives[, "name"]]/sols$exp_pen

  ## Drop any dominated solutions
  nondom_sols <- sols[apply(sols, 1, is_nondom, b=sols, objectives=objectives), ]
  ## check for duplicates
  sub <- unique(nondom_sols[, objectives[, "name"], drop=F])

  PS <- nondom_sols[row.names(sub),]
  PS
}
