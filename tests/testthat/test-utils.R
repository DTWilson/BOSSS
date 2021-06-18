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

test_that("expected improvement calculated", {
  set.seed(1)
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

  hypotheses <- data.frame(mu = c(0.3, 0),
                           var_u = c(0.05, 0.05),
                           var_e = c(0.95, 0.95))

  to_model <- data.frame(out_i = c(1),
                         hyp_i = c(1))

  get_det_obj <- function(design)
  {
    matrix(design, ncol = 2)[,1:2]
  }

  DoE <- init_DoE(20, design_space)

  DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=100, sim=sim_trial)))
  DoE$N <- 100

  models <- fit_models(DoE, to_model, design_space)

  PS <- best(design_space, models, DoE, objectives, constraints, to_model, get_det_obj)

  EI <- exp_improve(design = c(100, 10), N=100, PS, models,
                    design_space, constraints, objectives, get_det_obj, out_dim)

  expect_type(EI, "double")

  ptm <- proc.time()
  opt <- pso::psoptim(rep(NA, 2), exp_improve, lower=design_space$low, upper=design_space$up,
                      N=100, PS=PS, mod=models, design_space=design_space, constraints=constraints,
                      objectives=objectives, get_det_obj=get_det_obj, out_dim=3,
                      control=list(vectorize = T))
  proc.time() - ptm

  ptm <- proc.time()
  opt <- RcppDE::DEoptim(exp_improve, lower=design_space$low, upper=design_space$up,
                         control=list(trace=FALSE),
        N=100, PS=PS, mod=models, design_space=design_space, constraints=constraints,
        objectives=objectives, get_det_obj=get_det_obj, out_dim=3)
  proc.time() - ptm

  sol <- opt$par
  sol[1:2] <- round(sol[1:2])

  y <- calc_rates(sol, hypotheses=hypotheses, N=100, sim=sim_trial)

  DoE <- rbind(DoE, c(sol, y, 100))
})
