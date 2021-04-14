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

