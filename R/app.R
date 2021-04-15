library(shiny)

BOSSSapp <- function(...) {

  sim_trial <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]
    mu <- hypothesis[1]

    m <- n/k
    s_c <- sqrt(0.05 + 0.95/m)
    x0 <- rnorm(k, 0, s_c); x1 <- rnorm(k, mu, s_c)
    c(t.test(x0, x1)$p.value >= 0.05, n, k)
  }

  design_space <- data.frame(name = c("n", "k"),
                             low = c(100, 10),
                             up = c(500, 100),
                             int = c(TRUE, TRUE))

  hypotheses <- data.frame(matrix(c(0.3), nrow = 1))
  names(hypotheses) <- "mu"

  constraints <- data.frame(name = c("beta"),
                            out_i = c(1),
                            hyp_i = c(1),
                            nom = c(0.2),
                            delta = c(0.975),
                            stoch = c(TRUE)
  )

  objectives <- data.frame(name = c("f1", "f2"),
                           out_i = c(2, 3),
                           hyp_i = c(1, 1),
                           weight = c(2/5, 1),
                           stoch = c(FALSE, FALSE)
  )
  objectives$weight <- objectives$weight/sum(objectives$weight)
  objectives$name <- as.character(objectives$name)
  out_dim <- 3

  to_model <- data.frame(out_i = c(1),
                         hyp_i = c(1))

  ui <- fluidPage(
    # number of initial DoE
    numericInput("size", "Inital DoE size", value = 20),

    # number of MC evals
    numericInput("N", "Number of MC evals", value = 100),

    # Button to evaluate
    actionButton("initButton", "Initialise DoE"),

    tableOutput("table")
  )

  server <- function(input, output, session) {
    output$table <- renderTable({
      input$initButton

      DoE <- init_DoE(input$size, design_space)

      N <- 100
      DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=input$N, sim=sim_trial)))
      DoE$N <- input$N

      models <- fit_models(DoE, to_model, design_space)

      b <- best(design_space, models, DoE, objectives, constraints, to_model)
      b
    })
  }

  shinyApp(ui, server, ...)
}
