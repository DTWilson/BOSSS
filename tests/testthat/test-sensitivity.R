test_that("", {
  examples <- readRDS(test_path("examples", "examples.rds"))
  for(i in length(examples)){
    prob <- examples[[i]][[1]]

    p_name <- rownames(prob$hypotheses)[1]
    hyp_name <- colnames(prob$hypotheses)[1]
    p_val <- prob$hypotheses[p_name, hyp_name]

    sol <- examples[[i]][[2]]
    design <- sol$p_set[1,]

    m <- sensitivity(design, name = p_name, hypothesis = hyp_name,
                     lower = 0.5*p_val, upper = 1.5*p_val,
                     problem = prob, num_eval = 5, N = 20)

    expect_equal(nrow(m), 5)
  }
})
