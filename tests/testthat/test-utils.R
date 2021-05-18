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

  constraints <- data.frame(name = c("beta"),
                            out_i = c(1),
                            hyp_i = c(1),
                            nom = c(0.2),
                            delta = c(0.975),
                            stoch = c(T)
  )

  objectives <- data.frame(name = c("f1", "f2"),
                           out_i = c(2, 3),
                           hyp_i = c(1, 1),
                           weight = c(2/5, 1),
                           stoch = c(F, F)
  )
  objectives$weight <- objectives$weight/sum(objectives$weight)
  objectives$name <- as.character(objectives$name)
  out_dim <- 3

  to_model <- data.frame(out_i = c(1),
                         hyp_i = c(1))

  models <- fit_models(DoE, to_model, design_space)

  expect_equal(length(models), nrow(to_model))
  expect_is(models[[1]], "km")
})

test_that("mc estimates generated", {
  sim_trial <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]
    mu <- hypothesis[,1]; var_u <- hypothesis[,2]; var_e <- hypothesis[,3]

    m <- n/k
    s_c <- sqrt(var_u + var_e/m)

    x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
    c(stats::t.test(x0, x1)$p.value >= 0.05, n, k)
  }

  design <- c(400, 15)
  hypotheses <- t(c(0.3, 0.05, 0.95))

  results <- calc_rates(design, hypotheses, N = 50, sim = sim_trial)

  expect_equal(length(results), 2*3)
})
