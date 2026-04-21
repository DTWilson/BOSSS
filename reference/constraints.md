# Create a set of constraints

Create a set of constraints

## Usage

``` r
constraints(name, out, nom, delta, binary, hyp = NULL)
```

## Arguments

- name:

  character vector of constraint names.

- out:

  character vector denoting which simulation output each constraint
  pertains to.

- nom:

  numeric vector of nominal upper limits.

- delta:

  numeric vector of probabilities.

- binary:

  boolean vector denoting if the constraint output is binary or
  otherwise.

- hyp:

  optional character vector denoting which hypothesis each constraint
  pertains to. Defaults to NULL, in which case the simultion function's
  default hypothesis is used for all constraints.

## Value

A data.frame defining the constraints.

## Examples

``` r
constraints(name = c("tII"),
            out = c("s"),
            hyp = c("alt"),
            nom = c(0.1),
            delta = c(0.95),
            binary = c(TRUE))
#>   name out nom delta binary hyp
#> 1  tII   s 0.1  0.95   TRUE alt
```
