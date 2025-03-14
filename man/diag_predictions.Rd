% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{diag_predictions}
\alias{diag_predictions}
\title{Examine model predictions}
\usage{
diag_predictions(problem, solution, type = "response")
}
\arguments{
\item{problem}{BOSSS problem.}

\item{solution}{BOSSS solution.}

\item{type}{the type of prediction required for binary outcomes. The default
is on the scale of the response variable ("response"); the scale of the linear
predictor ("link") can be used instead.}
}
\value{
A list of dataframes, one for each model, giving the empirical
(Monte Carlo) point and interval estimates alongside their predicted point
and interval estimates, flagging when these do not agree.
}
\description{
Examine model predictions
}
