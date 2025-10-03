test_that("Solution can be initialised", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  N <- 50
  # Skipping example 3 (PRESSURE2) as too slow
  for(i in c(1,2,4)){
    prob <- examples[[i]][[1]]
    size <- 10*prob$dimen
    sol <- BOSSS_solution(size, N, prob)
    # Check Pareto set exisits
    expect_gt(nrow(sol$p_set), 0)
  }
})

test_that("Empirical estimates have correct dimensions", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]
    r <- PS_empirical_ests(sol, prob)
    expect_equal(ncol(r[[1]]), nrow(prob$objectives) + sum(prob$objectives$stoch))
    expect_equal(ncol(r[[2]]), nrow(prob$constraints) + sum(prob$constraints$stoch))
  }
})

