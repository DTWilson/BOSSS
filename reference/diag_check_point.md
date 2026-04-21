# Empirical check

Empirical check

## Usage

``` r
diag_check_point(design, problem, solution, N, current = NULL)
```

## Arguments

- design:

  Design to evaluate.

- problem:

  BOSSS problem.

- solution:

  BOSSS solution.

- N:

  Number of simulations to use when computing Monte Carlo estimates.

- current:

  An optional matrix containing the results of previous check_point
  calls, to be built upon.

## Value

A matrix with each row corresponds to a model of the BOSSS solution
object, giving the Monte Carlo estimates of the mean and variance along
side the number of simulations used to compute them.
