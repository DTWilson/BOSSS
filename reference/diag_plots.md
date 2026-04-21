# Plot one-dimensional model fits

Plot one-dimensional model fits

## Usage

``` r
diag_plots(design, problem, solution, type = "response")
```

## Arguments

- design:

  Design to centre plots at.

- problem:

  BOSSS problem.

- solution:

  BOSSS solution.

- type:

  the type of prediction required for binary outcomes. The default is on
  the scale of the response variable ("response"); the scale of the
  linear predictor ("link") can be used instead.

## Value

A list of plots of size (# models) x (# design variables).
