# Run a one-dimensional sensitivity analysis

Run a one-dimensional sensitivity analysis

## Usage

``` r
sensitivity(
  design,
  name = NULL,
  hypothesis = NULL,
  lower,
  upper,
  problem,
  num_eval = 20,
  N = 100
)
```

## Arguments

- design:

  Design to be tested

- name:

  name of the parameter to be varied

- hypothesis:

  name of the hypothesis around which the SA will be conducted

- lower:

  lower bound of the parameter

- upper:

  upper bound of the parameter

- problem:

  BOSSS problem

- num_eval:

  Number of points to evaluate

- N:

  Number of MC samples to use in each evaluation

## Value

A matrix of estimated means and their variances for each simulation
output over the sensitivity analysis parameter range.
