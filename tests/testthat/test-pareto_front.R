test_that("Pareto front dimensions are correct", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]

    pf <- pareto_front(sol, prob)[[1]]
    expect_equal(ncol(pf) - 1, nrow(prob$objectives))
  }
})

test_that("Constraint penalties are correct", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]
    for(j in 1:nrow(prob$constraints)){
      cons <- check_constraint(j, sol, prob)
      expect_true(all((cons[,2] >= prob$constraints$delta[j]) == (cons[,1] > 0.5)))
    }
  }
})



