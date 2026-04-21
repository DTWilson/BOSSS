# Examine model predictions

Examine model predictions

## Usage

``` r
diag_predictions(problem, solution, type = "response")
```

## Arguments

- problem:

  BOSSS problem.

- solution:

  BOSSS solution.

- type:

  the type of prediction required for binary outcomes. The default is on
  the scale of the response variable ("response"); the scale of the
  linear predictor ("link") can be used instead.

## Value

A list of dataframes, one for each model, giving the empirical (Monte
Carlo) point and interval estimates alongside their predicted point and
interval estimates, flagging when these do not agree.
