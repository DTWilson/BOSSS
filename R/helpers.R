# Some helper functions to create the four data frames which partly specify
# the problem.

design_space <- function(name, lower, upper) {
  stopifnot(upper > lower)

  data.frame(name = name,
             lower = lower,
             upper = upper)
}

hypotheses <- function(name, param_matrix) {
  df <- as.data.frame(param_matrix)
  row.names(df) <- name
  names(df) <- paste0("hyp_", 1:ncol(param_matrix))
  df
}

constraints <- function(name, out_i, hyp_i, nom, delta, stoch) {

  data.frame(name = name,
             out_i= out_i,
             hyp_i = hyp_i,
             nom = nom,
             delta = delta,
             stoch = stoch)
}

objectives <- function(name, out_i, hyp_i, weight, stoch) {

  data.frame(name = name,
             out_i= out_i,
             hyp_i = hyp_i,
             weight = weight,
             stoch = stoch)
}
