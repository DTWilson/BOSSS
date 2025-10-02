test_that("Solution can be initialised", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  N <- 10
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    size <- 10*prob$dimen
    sol <- BOSSS_solution(size, N, prob)
    # Check Pareto set exisits
    expect_gt(nrow(sol$p_set), 0)
  }
})
