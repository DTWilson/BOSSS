test_that("initial DoE created", {
  design_space <- data.frame(name = c("n", "k"),
                             low = c(100, 10),
                             up = c(500, 100),
                             int = c(TRUE, TRUE))
  size <- 10

  DoE <- init_DoE(size, design_space)

  expect_equal(size, nrow(DoE))
})


test_that("models fit", {
  design_space <- data.frame(name = c("n", "k"),
                             low = c(100, 10),
                             up = c(500, 100),
                             int = c(TRUE, TRUE))

  DoE <- data.frame(n = c(300, 400, 200),
                    k = c(55, 32, 78),
                    a = c(0.11, 0.10, 0.15),
                    b = c(0.000988, 0.000909, 0.00128),
                    N = c(100, 100, 100))

  to_model <- data.frame(out_i = c(1),
                         hyp_i = c(1))

  models <- fit_models(DoE, to_model, design_space)

  expect_equal(length(models), nrow(to_model))
  expect_is(models[[1]], "km")
})

