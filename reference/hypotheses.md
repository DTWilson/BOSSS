# Create a set of hypotheses

Create a set of hypotheses

## Usage

``` r
hypotheses(par_name = NULL, sim = NULL, values, hyp_names)
```

## Arguments

- par_name:

  optional character vector of model parameter names.

- sim:

  optional simulation function.

- values:

  numeric matrix, each column giving the model parameter values under a
  specific hypothesis.

- hyp_names:

  character vector of hypothesis names.

## Value

A data.frame defining the hypotheses.

## Examples

``` r
hypotheses(values = matrix(c(0.3, 0.05, 0.95), ncol = 1),
           hyp_names = c("alt"),
           par_name = c("mu", "var_u", "var_e"))
#>        alt
#> mu    0.30
#> var_u 0.05
#> var_e 0.95
```
