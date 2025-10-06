test_that("problems can be constructed", {

  sim_cRCT <- function(design = list(m = 10, n = 20),
                       hypothesis = list(beta_1 = 0.3, var_e = 0.95, var_u = 0.05)){

    m <- design[[1]]; n <- design[[2]]
    beta_1 <- hypothesis[[1]]; var_e <- hypothesis[[2]]; var_u <- hypothesis[[3]]

    m <- floor(m); n <- floor(n)

    x <- rep(c(0,1), each = m)
    y <- rnorm(2*m, 0, sqrt(var_u + var_e/n))
    y <- y + x*beta_1

    s <- stats::t.test(y[x==1], y[x==0], alternative = "greater")$p.value > 0.025

    return(c(s = s))
  }

  det_cRCT <- function(design = list(m = 10, n = 20),
                       hypothesis = list(beta_1 = 0.3, var_e = 0.95, var_u = 0.05))
  {
    m <- design[[1]]; n <- design[[2]]
    beta_1 <- hypothesis[[1]]; var_e <- hypothesis[[2]]; var_u <- hypothesis[[3]]

    return(c(m = floor(m), N = floor(n)*floor(m)))
  }

  design_space <- design_space(lower = c(10, 5),
                               upper = c(50, 20),
                               sim = sim_cRCT)

  hypotheses <- hypotheses(values = matrix(c(0.3, 0.95, 0.05), ncol = 1),
                           hyp_names = c("alt"),
                           sim = sim_cRCT)

  constraints <- constraints(name = c("con_tII"),
                             out = c("s"),
                             hyp = c("alt"),
                             nom = c(0.2),
                             delta = c(0.95),
                             binary = c(TRUE))

  objectives <- objectives(name = c("N", "m"),
                           out = c("N", "m"),
                           hyp = c("alt", "alt"),
                           weight = c(1, 10),
                           binary = c(FALSE, FALSE))

  prob <- BOSSS_problem(sim_cRCT, design_space, objectives, hypotheses, constraints, det_func = det_cRCT)

  expect_true(is(prob, "BOSSS_problem"))
})
