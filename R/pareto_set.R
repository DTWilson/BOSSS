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
pareto_set <- function(design_space, models, DoE, objectives, constraints, to_model, get_det_obj, b=NULL)
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
