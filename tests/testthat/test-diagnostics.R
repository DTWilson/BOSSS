test_that("Correct number of plots are returned", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]
    design <- sol$p_set[1,]
    plots <- diag_plots(design, prob, sol, type = "response")
    num_models <- nrow(sol$to_model)
    expect_equal(length(plots), num_models*prob$dimen)

    plots <- diag_plots(design, prob, sol, type = "link")
    num_models <- nrow(sol$to_model)
    expect_equal(length(plots), num_models*prob$dimen)
  }
})

test_that("Empirical intervals returned", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]
    design <- sol$p_set[1,]
    num_models <- nrow(sol$to_model)

    r <- diag_check_point(design, prob, sol, N=50)
    expect_length(r, 3*num_models)
  }
})

test_that("Prediction diags returns a dataframe of correct size", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]

    dfs <- diag_predictions(prob, sol, type = "response")
    for(j in 1:length(dfs)){
      expect_equal(ncol(dfs[[j]]), prob$dimen + 7)
    }
  }
})
