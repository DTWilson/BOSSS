BOSSSapp <- function(...) {

  sim_trial <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]
    mu <- hypothesis[,1]; var_u <- hypothesis[,2]; var_e <- hypothesis[,3]

    m <- n/k
    s_c <- sqrt(var_u + var_e/m)
    x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
    c(stats::t.test(x0, x1)$p.value >= 0.05, n, k)
  }

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

  to_model <- data.frame(out_i = c(1),
                         hyp_i = c(1))

  ui <- shiny::fluidPage(

    # Design space matrix
    shinyMatrix::matrixInput("DSnums", class = "numeric",
                cols = list(names = TRUE), rows = list(extend = TRUE, names = TRUE, editableNames = TRUE),
                value =  matrix(c(100, 10, 500, 100, 1, 1), 2, 3,
                                dimnames = list(c("n", "k"), c("Min", "Max", "Integer")))),

    # Hypotheses matrix
    shinyMatrix::matrixInput("Hypnums", class = "numeric",
                             cols = list(extend = TRUE, names = TRUE, editableNames = TRUE),
                             rows = list(extend = TRUE, names = TRUE, editableNames = TRUE),
                             value =  matrix(c(0.3, 0, 0.05, 0.05, 0.95, 0.95), ncol = 3,
                                             dimnames = list(letters[1:2], c("mu", "var_u", "var_e")))),

    # number of initial DoE
    shiny::numericInput("size", "Inital DoE size", value = 20),

    # number of MC evals
    shiny::numericInput("N", "Number of MC evals", value = 100),

    # Button to evaluate
    shiny::actionButton("initButton", "Initialise DoE"),

    shiny::tableOutput("table"),

    shiny::verbatimTextOutput("code")
  )

  server <- function(input, output, session) {

    ds <- shiny::reactive(
      data.frame(name = rownames(input$DSnums),
                              low = input$DSnums[,1],
                              up = input$DSnums[,2],
                              int = input$DSnums[,3])[1:(nrow(input$DSnums) - 1),]
    )

    hyps <- shiny::reactive(
      data.frame(input$Hypnums)
      )

    output$table <- shiny::renderTable({
      input$initButton

      design_space <- ds()

      m <- hyps()
      m <- m[rowSums(is.na(m)) != ncol(m),]
      m <- m[, colSums(is.na(m)) != nrow(m)]
      #hypotheses <- t(c(0.3, 0.05, 0.95))

      DoE <- init_DoE(input$size, design_space)

      N <- 100
      DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=input$N, sim=sim_trial)))
      DoE$N <- input$N
      DoE

      models <- fit_models(DoE, to_model, design_space)

      b <- best(design_space, models, DoE, objectives, constraints, to_model)
      b
    })

    output$code <- shiny::renderPrint({
      ds()
      #hyps()
    })
  }

  shiny::shinyApp(ui, server, ...)
}
