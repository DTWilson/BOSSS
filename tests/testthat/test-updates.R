test_that("can update constraints", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]

    prob <- update_constraint(prob, number = 1, nom = 2*prob$constraints$nom[1])
    sol <- update_solution(sol, prob)

    # Constraint relaxed so would expect more valid designs
    expect_gte(nrow(sol$p_set), 1)
  }
})

test_that("can extend initial solution", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]

    old_size <- nrow(sol$DoE)

    sol <- extend_initial(prob, sol, extra_N = 10, extra_points = 1)

    expect_equal(nrow(sol$DoE), old_size + 1)
  }
})
