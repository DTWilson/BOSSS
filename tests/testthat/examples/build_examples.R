# Take example problems and solutions from our vignettes and save them
# together to be used when testing

examples <- list()

prob <- readRDS(here("vignettes", "example_files", "cRCT_prob.rds"))
sol <- readRDS(here("vignettes", "example_files", "cRCT_sol_final.rds"))

examples[[1]] <- list(prob, sol)

prob <- readRDS(here("vignettes", "example_files", "NIFTy_prob.rds"))
sol <- readRDS(here("vignettes", "example_files", "NIFTy_sol_final.rds"))

examples[[2]] <- list(prob, sol)

prob <- readRDS(here("vignettes", "example_files", "P2_prob.rds"))
sol <- readRDS(here("vignettes", "example_files", "P2_sol_final.rds"))

examples[[3]] <- list(prob, sol)

saveRDS(examples, test_path("examples", "examples.rds"))


