test_that("Initial DoE created", {
  design_space <- data.frame(name = c("n", "k"),
                             low = c(100, 10),
                             up = c(500, 100),
                             int = c(TRUE, TRUE))
  size <- 10

  DoE <- init_DoE(size, design_space)

  expect_equal(size, nrow(DoE))
})
