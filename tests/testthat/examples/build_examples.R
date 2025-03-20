# Setup some example BOSSS_problems and BOSSS_simulations which we can and then
# use in our unit tests

examples <- list()

###########################################
## 1. Cluster RCT from the main vignette ##
###########################################

sim_cRCT <- function(k = 10, n = 200, beta_1 = 0.3, sig_w = 0.975, sig_b = 0.224){
  k <- floor(k); n <- floor(n)
  u_c <- rnorm(k, sd = sig_b)
  cluster_size <-  rep(n%/%k, k) + c(rep(1, n%%k), rep(0, k - n%%k))
  y_c <- rep(u_c, times = cluster_size) + rnorm(n, 0, sig_w)

  u_i <- rnorm(k, sd = sig_b)
  y_i <- beta_1 + rep(u_i, times = cluster_size) + rnorm(n, 0, sig_w)

  s <- t.test(y_i, y_c, alternative = "greater")$p.value > 0.05

  return(c(s = s))
}

det_cRCT <- function(k = 10, n = 200, beta_1 = 0.3, sig_w = 0.975, sig_b = 0.224){
  return(c(k = k, n = n, m = n/k))
}

design_space <- design_space(lower = c(10, 100),
                             upper = c(100, 500),
                             sim = sim_cRCT)

hypotheses <- hypotheses(values = matrix(c(0.3, sqrt(0.95), sqrt(0.05)), ncol = 1),
                         hyp_names = c("alt"),
                         sim = sim_cRCT)

constraints <- constraints(name = c("con_tII", "con_m"),
                           out = c("s", "m"),
                           hyp = c("alt", "alt"),
                           nom = c(0.2, 10),
                           delta = c(0.95, 1),
                           binary = c(TRUE, FALSE))

objectives <- objectives(name = c("n", "k"),
                         out = c("n", "k"),
                         hyp = c("alt", "alt"),
                         weight = c(1, 10),
                         binary = c(FALSE, FALSE))

prob <- BOSSS_problem(sim_cRCT, design_space, hypotheses, objectives, constraints, det_func = det_cRCT)

size <- 20
N <- 500
sol <- BOSSS_solution(size, N, prob)

examples[[1]] <- list(prob, sol)

###############################################
## 2. Two stage trial from examples vignette ##
###############################################

sim_2S <- function(n1 = 100, c1 = 0.1, n2 = 100, c2 = 0.2,
                   delta = 0.3, sd = 1){
  y0 <- rnorm(n1, 0, sd)
  y1 <- rnorm(n1, delta, sd)

  go <- (mean(y1) - mean(y0)) > c1

  if(go){
    y0 <- c(y0, rnorm(n2, 0, sd))
    y1 <- c(y1, rnorm(n2, delta, sd))

    go2 <- (mean(y1) - mean(y0)) > c2
    if(go2){
      g <- TRUE
      s <- FALSE
    } else {
      g <- FALSE
      s <- TRUE
    }
    n <- 2*(n1 + n2)

  } else {
    g <- FALSE
    s <- TRUE
    n <- 2*n1
  }

  return(c(g = g, s = s, n = n))
}

design_space <- design_space(lower = c(10, -5, 10, -5),
                             upper = c(500, 5, 500, 5),
                             sim = sim_2S)

hypotheses <- hypotheses(values = matrix(c(0, 1, 0.3, 1), ncol = 2),
                         hyp_names = c("null", "alt"),
                         sim = sim_2S)

constraints <- constraints(name = c("con_tI", "con_tII"),
                           out = c("g", "s"),
                           hyp = c("null", "alt"),
                           nom = c(0.2, 0.4),
                           delta = c(0.95, 0.95),
                           binary = c(TRUE, TRUE))

objectives <- objectives(name = c("tI", "tII", "En"),
                         out = c("g", "s", "n"),
                         hyp = c("null", "alt", "null"),
                         weight = c(100, 100, 1),
                         binary = c(TRUE, TRUE, FALSE))

prob <- BOSSS_problem(sim_2S, design_space, hypotheses, objectives, constraints)

size <- 40
N <- 200
sol <- BOSSS_solution(size, N, prob)

for(i in 1:10) {
  sol <- iterate(sol, prob, N)
}

examples[[2]] <- list(prob, sol)

#saveRDS(examples, test_path("examples", "examples.rds"))


