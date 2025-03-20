test_that("itaration updates solution", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]

    new_sol <- iterate(sol, prob, N = 100)

    expect_equal(nrow(new_sol$results[[1]]), nrow(sol$results[[1]]) + 1)
  }
})
