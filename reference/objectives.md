# Create a set of objectives

Create a set of objectives

## Usage

``` r
objectives(name, out, weight, binary = NULL, hyp = NULL)
```

## Arguments

- name:

  character vector of objective names.

- out:

  character vector denoting which simulation output each objective
  pertains to.

- weight:

  numeric vector of weights assigned to each objective.

- binary:

  optional boolean vector denoting if the output of the objective
  function is binary (TRUE) or continuous (FALSE).

- hyp:

  optional character vector denoting which hypothesis each constraint
  pertains to. Defaults to NULL, in which case the simultion function's
  default hypothesis is used for all objectives.

## Value

A data.frame defining the objectives.

## Examples

``` r
objectives(name = c("min_n", "min_k"),
                     out = c("n", "k"),
                     hyp = c("alt", "alt"),
                     weight = c(10, 1))
#>    name out weight hyp
#> 1 min_n   n     10 alt
#> 2 min_k   k      1 alt
```
