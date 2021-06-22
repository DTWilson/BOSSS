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
#' hypotheses <- t(c(0.3))
#' calc_rates(design, hypotheses, N = 100, sim = sim_trial)
calc_rates <- function(design, hypotheses, N, sim)
{
  results <- NULL
  for(i in 1:nrow(hypotheses)){
    sims <- replicate(N, sim(design, as.data.frame(hypotheses)[i,]))
    for(j in 1:nrow(sims)){
      results <- c(results, mean(sims[j,]), stats::var(sims[j,])/N)
      names(results)[(length(results) -1)] <- paste0(rownames(sims)[j], "_m_", rownames(hypotheses)[i])
      names(results)[length(results)] <- paste0(rownames(sims)[j], "_v_", rownames(hypotheses)[i])
    }
  }
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
best <- function(design_space, models, DoE, objectives, constraints, to_model, get_det_obj, b=NULL)
{
  ## Return the set of current Pareto optimal solutions,
  ## penalising constrain violations and considering only solutions
  ## where some evaluation has actually happened
  sols <- DoE

  dim <- nrow(design_space)
  out_dim <- max(c(objectives$out_i, constraints$out_i))

  ## Get objective values
  obj_v <- t(apply(sols[,1:2], 1, function(x) predict_obj(x, models, objectives, get_det_obj, dim, to_model)))
  sols <- cbind(sols, f1=obj_v[,1], f2=obj_v[,2])

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

predict_obj <- function(design, models, objectives, get_det_obj, dim, to_model)
{
  obj_vals <- NULL
  i <- 1
  for(i in 1:nrow(objectives)){
    if(objectives$stoch[i]){
      ## For stochastic objectives, we are just taking the lowest quantile and
      ## treating this as deterministic - need to improve
      model_index <- which(to_model$out_i == objectives[i, "out_i"] &
                             to_model$hyp_i == objectives[i, "hyp_i"])
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design[,1:dim, drop=F], type="SK")
      f <- p$mean - qnorm(0.7)*p$sd
      obj_vals <- c(obj_vals, f*objectives$weight)
    } else {
      ## For deterministic objectives, use the user-written function
      ## For now, assume independent of hypothesis
      f <- get_det_obj(design)[objectives[i, "out_i"]]
      obj_vals <- c(obj_vals, f*objectives$weight[i])
    }
  }
  t(objectives$weight*t(matrix(obj_vals, ncol = nrow(objectives))))
}

exp_improve <- function(design, N, PS, models, design_space, constraints, objectives, get_det_obj, out_dim, to_model)
{
  dim <- nrow(design_space)

  ## Get objective value of design
  fs <- predict_obj(design, models, objectives, get_det_obj, dim)

  design <- as.data.frame(matrix(design, ncol = dim))
  names(design) <- design_space$name

  ## Get expected penalisation if we were to evaluate at design,
  ## using the models of constraint functions
  design$exp_pen <- 1
  for(i in 1:nrow(constraints)){
    index <- DoE_index(constraints[i, "out_i"], constraints[i, "hyp_i"],
                       dim, out_dim)
    model_index <- which(to_model$out_i == constraints[i, "out_i"] & to_model$hyp_i == constraints[i, "hyp_i"])
    nom <- constraints[i, "nom"]
    if(constraints[i, "stoch"]){
      p <- DiceKriging::predict.km(models[[model_index]], newdata=design[,1:dim, drop=F], type="SK")
      mc_vars <- 0.25/N #p$mean*(1-p$mean)/N
      pred_q_mean <- p$mean + stats::qnorm(constraints[i, "delta"])*sqrt(mc_vars*(p$sd^2)/(mc_vars+(p$sd^2)))
      pred_q_var <- ((p$sd^2)^2)/(mc_vars+(p$sd^2))
      design$exp_pen <- design$exp_pen*stats::pnorm(nom, pred_q_mean, sqrt(pred_q_var))
    } else {
      pen <- ifelse(DoE[, index] > nom, 0.0000001, 1)
    }
  }

  ## Improvement is quantified by the increase in the dominated hypervolume of
  ## the approximation set if this design was included
  ref <- design_space$up*objectives$weight
  PS2 <- as.matrix(PS[, objectives$name])
  current <- mco::dominatedHypervolume(PS2, ref)
  pos <- apply(fs, 1, function(obj) mco::dominatedHypervolume(as.matrix(rbind(PS2, obj)), ref) )
  #pos <- mco::dominatedHypervolume(as.matrix(rbind(PS2, fs)), ref)
  imp <- (current-pos)*design$exp_pen

  ## Minimising, so keeping negative
  return(imp)
}
