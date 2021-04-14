#' Bind factors
#'
#' Some description
#'
#' @param a factor
#' @param b factor
#'
#' @return factor
#' @export
#'
#' @examples
#' fbind(iris$Species[c(1, 51, 101)], PlantGrowth$group[c(1, 11, 21)])
fbind <- function(a, b) {
  factor(c(as.character(a), as.character(b)))
}

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

