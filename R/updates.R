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
