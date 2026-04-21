# Perform one iteration of Bayesian optimisation

Perform one iteration of Bayesian optimisation

## Usage

``` r
iterate(solution, problem, N, design = NULL)
```

## Arguments

- solution:

  current BOSSS solution.

- problem:

  BOSSS problem.

- N:

  number of simulations to use when computing Monte Carlo estimates.

- design:

  optional vector in the design space to be evaluated. If left NULL, an
  optimal design will be sought.

## Value

An updated BOSSS solution object.
