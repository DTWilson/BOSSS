test_that("DoE initialisation works", {
  size <- 10
  ds <- design_space(name = c("a", "b"),
                     lower = c(-10, 30),
                     upper = c(50, 40))
  expect_equal(nrow(init_DoE(size, ds)), size)
})
