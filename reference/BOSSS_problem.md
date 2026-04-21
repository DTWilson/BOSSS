# Create a BOSSS problem

Create a BOSSS problem

## Usage

``` r
BOSSS_problem(
  sim_trial,
  design_space,
  objectives,
  hypotheses = NULL,
  constraints = NULL,
  det_func = NULL
)
```

## Arguments

- sim_trial:

  function which generates a single (possibly multivariate) Monte Carlo
  outcome of a design under a hypothesis.

- design_space:

  data frame constructed via
  [`design_space()`](https://dtwilson.github.io/BOSSS/reference/design_space.md).

- objectives:

  data frame constructed via
  [`objectives()`](https://dtwilson.github.io/BOSSS/reference/objectives.md).

- hypotheses:

  data frame constructed via
  [`hypotheses()`](https://dtwilson.github.io/BOSSS/reference/hypotheses.md).

- constraints:

  optional ata frame constructed via
  [`constraints()`](https://dtwilson.github.io/BOSSS/reference/constraints.md).

- det_func:

  optional function which generates deterministic outcomes of a design
  under a hypothesis.

## Value

An object of class BOSSS_problem.
