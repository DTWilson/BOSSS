BOSSSapp <- function(...) {

  sim_trial <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]
    mu <- hypothesis[1]

    m <- n/k
    s_c <- sqrt(0.05 + 0.95/m)
    x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
    c(stats::t.test(x0, x1)$p.value >= 0.05, n, k)
  }

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

  to_model <- data.frame(out_i = c(1),
                         hyp_i = c(1))

  #get_design_var <- function(id)
  #{
  #  list(textInput(paste0(id, "Lab"), label = "Design variable name", value = id),
  #       numericInput(paste0(id, "Min"), "Minimum", value = 10),
  #       numericInput(paste0(id, "Max"), "Maximum", value = 100),
  #       numericInput(paste0(id, "Int"), "Integer", value = 1))
  #}

  ui <- shiny::fluidPage(

    shinyMatrix::matrixInput("DSnums", class = "numeric",
                cols = list(names = TRUE), rows = list(extend = TRUE, names = TRUE, editableNames = TRUE),
                value =  matrix(c(100, 10, 500, 100, 1, 1), 2, 3,
                                dimnames = list(c("n", "k"), c("Min", "Max", "Integer")))),

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

    #ds <- reactive(data.frame(name = c(input$aLab, input$bLab),
    #                                    low = c(input$aMin, input$bMin),
    #                                    up = c(input$aMax, input$bMax),
    #                                    int = c(input$aInt==1, input$bInt==1)))

    ds <- shiny::reactive(
      data.frame(name = rownames(input$DSnums), #unlist(strsplit(input$DSnames,",")),
                              low = input$DSnums[,1],
                              up = input$DSnums[,2],
                              int = input$DSnums[,3])[1:(nrow(input$DSnums) - 1),]
      )

    output$table <- shiny::renderTable({
      input$initButton

      design_space <- ds()

      DoE <- init_DoE(input$size, design_space)

      N <- 100
      DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=input$N, sim=sim_trial)))
      DoE$N <- input$N

      models <- fit_models(DoE, to_model, design_space)

      b <- best(design_space, models, DoE, objectives, constraints, to_model)
      b
    })

    output$code <- shiny::renderPrint({
      ds()
    })
  }

  shiny::shinyApp(ui, server, ...)
}
