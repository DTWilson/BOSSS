test_that("Penalties are bounded", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in 1:length(examples)){
    prob <- examples[[i]][[1]]
    sol <- examples[[i]][[2]]
    design <- sol$p_set[1,][1,1:prob$dimen]

    pen <- exp_penalty(design, prob, sol, N=100)

    expect_gte(pen, 0)
    expect_lte(pen, 1)
  }
})
