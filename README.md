
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BOSSS <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/DTWilson/BOSSS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/DTWilson/BOSSS/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of BOSSS is to help people use Bayesain optimisation to solve
sample size determination problems when simulation is required to
calculate operating characteristics.

## Installation

You can install the development version of BOSSS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DTWilson/BOSSS")
```

## Example

Suppose we want to design a cluster randomised RCT which will compare
the means of our two groups via a t test. Specifically, we want to find
the smallest trial which have a power of 90%, where power must be
estimated numerically using simulation. To solve this problem using
BOSSS we start by writing two functions:

``` r
library(BOSSS)

sim_trial <- function(n=200, k=10, mu=0.3, var_u=0.05, var_e=0.95){

  # Number of patients in each cluster
  m <- n/k
  
  # SD of cluster means
  s_c <- sqrt(var_u + var_e/m)
  
  # Simulate cluster means
  x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
  
  return(c(s = stats::t.test(x0, x1)$p.value >= 0.05))
}

det_func <- function(n=200, k=10, mu=0.3, var_u=0.05, var_e=0.95){
  return(c(n=n, k=k))
}
```

The first function is a simulation, returning a binary output indicating
if the t test failed to return a significant result. The second function
is deterministic, returning the per-arm sample size and per-arm number
of clusters. There are some conventions which these functions must
adhere to:

i The argument lists of the two functions must be identical; ii
Arguments should start with all design variables, followed by all model
parameters; iii Design variables must be continuous; iv Outputs must be
a named vectors; v The simulation function is mandatory, but the
determinsitc function is optional.

In this example, the design variables are the total per-arm sample size
`n` and the number of clusters in each arm, `k`. We put ranges on these
to define our **design space**:

``` r
design_space <- design_space(lower = c(10, 3),
                             upper = c(500, 50),
                             sim = sim_trial)
```

We have decided to search over a space of trial designs with a total
per-arm sample size between 10 and 500, and a per-arm number of clusters
from 3 to 50. Note that we are extracting the names of our design
variables from the simulation function. Alternatively, we can specify
them manually using the `names` argument.

The model parameters are the mean difference `mu`, the between-cluster
variance `var_u`, and the within-cluster variance `var_e`. We choose
which parameter values we will be simulating under and call these sets
**hypotheses**:

``` r
hypotheses <- hypotheses(values = matrix(c(0.3, 0.05, 0.95), ncol = 1),
                         hyp_names = c("alt"),
                         sim = sim_trial)
```

To estimate power we will simulate under a single alternative
hypothesis. We want to impose a nominal level on this; in particular, we
want the type II error rate to be less than 0.1. We specify this as a
**constraint**:

``` r
constraints <- constraints(name = c("tII"),
                   out = c("s"),
                   hyp = c("alt"),
                   nom = c(0.1),
                   delta = c(0.95),
                   binary = c(TRUE))
```

Each constraint must be tied to a hypothesis and an output. In this
case, the output is `s`, which was defined in the simulation function as
a binary event indicating if the trial failed to return a statistically
significant result. The `delta` part of a constraint is the tolerance we
allow for it; in this case, we consider the constraint satisfied if
there is at least 0.95 chance that the type II error rate is less than
0.1.

Finally, we formalise the **objectives** we want to minimise:

``` r
objectives <- objectives(name = c("f1", "f2"),
                 out = c("n", "k"),
                 hyp = c("alt", "alt"),
                 weight = c(10, 1),
                 binary = c(FALSE, FALSE))
```

As with constraints, objectives are defined with respect to an output
and a hypothesis. In this case, we wish to minimise both the total
sample size and the number of clusters, both of which are outputs of the
deterministic function. Although BOSSS uses a non-scalarising approach
to solving multi-objective problems, we nevertheless require a relative
weighting for use in the internal optimisation process. Here, we specify
that minimising the number of clusters is around 10 times as valuable as
minimising the number of participant.

Together, these are the ingredients of a **problem** object:

``` r
problem <- BOSSS_problem(sim_trial, design_space, hypotheses, objectives, constraints, det_func = det_func)
```

With our problem specified, we generate a **solution** by evaluating a
first set of 20 possible designs, evenly spread over our design space,
using `N` simulations each time. Printing the result will show us the
estimated Pareto set (that is, the set of non-dominated designs):

``` r
size <- 20
N <- 500

solution <- BOSSS_solution(size, N, problem)
#> Checking simulation speed...
#> Initialisation will take approximately 1.365187 secs 
#> Solution found
print(solution)
#> Design variables for the Pareto set: 
#> 
#>          n       k
#> 9  346.875 41.1875
#> 13 408.125 35.3125
#> 
#> Corresponding objective function values... 
#> 
#>    n, alt (mean) k, alt (mean)
#> 9        346.875       41.1875
#> 13       408.125       35.3125
#> 
#> ...and constraint function values:
#> 
#>    s, alt (mean) s, alt (var)
#> 9      -2.493162   0.02834157
#> 13     -2.895870   0.04027675
```

To improve the solution we can **iterate** as many times as we like,
where each iteration will try to select the design which will give us
the biggest improvement and evaluate it.

``` r
for(i in 1:10){
  solution <- iterate(solution, problem, N)
}

print(solution)
#> Design variables for the Pareto set: 
#> 
#>           n        k
#> 9  346.8750 41.18750
#> 22 453.4098 29.56320
#> 23 499.9935 26.08903
#> 28 331.7451 43.80545
#> 29 387.7225 33.26920
#> 
#> Corresponding objective function values... 
#> 
#>    n, alt (mean) k, alt (mean)
#> 9       346.8750      41.18750
#> 22      453.4098      29.56320
#> 23      499.9935      26.08903
#> 28      331.7451      43.80545
#> 29      387.7225      33.26920
#> 
#> ...and constraint function values:
#> 
#>    s, alt (mean) s, alt (var)
#> 9      -2.493162   0.02834157
#> 22     -2.644208   0.03226089
#> 23     -2.465155   0.02767846
#> 28     -2.309640   0.02432073
#> 29     -2.261940   0.02339297
```

We can also visualise our solution by plotting the Pareto front (that
is, the objective values of the solutions in the Pareto set):

``` r
plot(solution)
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

After we have finished iterating, we can then select a specific design
from the estimated Pareto set. We may wish to double check the operating
characteristics of the design using a large number of simulations. We
can do this via BOSSS, which will print both the empirical 95%
confidence interval and the corresponding interval given by the Gaussian
Process surrogate model as a diagnostic.

``` r
design <- solution$p_set[1,]

r <- diag_check_point(design, problem, solution, N=10^5)
#> Model 1 prediction interval: [0.074, 0.095]
#> Model 1 empirical interval: [0.084, 0.087]
```
