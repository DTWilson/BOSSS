#' Get the Pareto front
#'
#' @param design_space data.frame
#' @param models list of objects of class km
#' @param DoE data.frame
#' @param objectives data.frame
#' @param constraints data.frame
#' @param to_model data.frame
#' @param out_dim
#' @param det_obj
#' @param b data.frame
#'
#' @return data.frame
#' @export
pareto_front <- function(models, DoE, results, to_model, problem, b=NULL)
{
  ## Return the objective values of current Pareto optimal solutions,
  ## penalising constrain violations and considering only solutions
  ## where some evaluation has actually happened
  dimen <- problem$dimen

  ## Get objective values
  # Add the objective values for all evaluated points
  obj_v <- predict_obj(DoE[,1:dimen], models, problem$objectives, problem$det_obj, to_model)

  ## Penalise constraint violations
  ## Note - some repetition with penalty in EHI so try to consolidate
  constraints <- problem$constraints
  exp_pen <- 1
  for(i in 1:nrow(constraints)){
    #index <- DoE_index(constraints[i, "out_i"], constraints[i, "hyp_i"], dim, out_dim)
    model_index <- which(to_model$out_i == constraints[i, "out_i"] & to_model$hyp_i == constraints[i, "hyp_i"])
    nom <- constraints[i, "nom"]
    if(constraints[i, "stoch"]){
      p <- DiceKriging::predict.km(models[[model_index]], newdata=DoE[,1:dimen, drop=F], type="UK")
      pen <- stats::pnorm(nom, p$mean, p$sd)
      pen <- ifelse(pen < constraints[i, "delta"], 0.0000001, 1)
      exp_pen <- exp_pen*pen
    } else {
      pen <- ifelse(results[[constraints[i, "hyp_i"], constraints[i, "out_i"]]][,1] > nom, 0.0000001, 1)
      exp_pen <- exp_pen*pen
    }
  }
  obj_v <- obj_v/exp_pen + 1/exp_pen - 1

  # add an ID variable for extracting the PS later
  obj_v <- cbind(obj_v, sapply(1:nrow(obj_v), function(x) which(order(obj_v[, ncol(obj_v)]) == x)))

  ## Get the Pareto front
  pf <- t(emoa::nondominated_points(t(obj_v)))
  return(list(unique(pf), obj_v[,ncol(obj_v)]))
}

#m <- t(nondominated_points(t(as.matrix(sols[,c("f1", "f2")]))))
#z <- sols[match(m[,1], sols[,"f1"]),]

