% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{init_DoE}
\alias{init_DoE}
\title{Set up initial space-filling design}
\usage{
init_DoE(size, design_space)
}
\arguments{
\item{size}{integer}

\item{design_space}{data.frame}
}
\value{
data.frame
}
\description{
Set up initial space-filling design
}
\examples{
design_space <- data.frame(name = c("n", "k"),
low = c(100, 10),
up = c(500, 100),
int = c(TRUE, TRUE))

init_DoE(10, design_space)
}
