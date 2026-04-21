# Create a design space

Create a design space

## Usage

``` r
design_space(name = NULL, sim = NULL, lower, upper)
```

## Arguments

- name:

  optional character vector of design variable names.

- sim:

  optional simulation function.

- lower:

  numeric vector of lower limits.

- upper:

  numeric vector of upper limits.

## Value

A data.frame defining the design space.

## Examples

``` r
design_space(lower = c(10, 3),
             upper = c(500, 50),
             name = c("n", "k"))
#>   name lower upper
#> 1    n    10   500
#> 2    k     3    50
```
