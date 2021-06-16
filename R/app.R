BOSSSapp <- function(...) {

  sim_trial <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]
    mu <- hypothesis[,1]; var_u <- hypothesis[,2]; var_e <- hypothesis[,3]

    m <- n/k
    s_c <- sqrt(var_u + var_e/m)
    x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
    c(s = stats::t.test(x0, x1)$p.value >= 0.05, p = n, c = k)
  }

  det_out <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]

    c(s = NA, p = n, c = k)
  }

  #num_hyp <- 2

  ui <- shiny::fluidPage(

    # Design space matrix
    shinyMatrix::matrixInput("DSnums", class = "numeric",
                cols = list(names = TRUE), rows = list(names = TRUE),
                value =  matrix(c(100, 10, 500, 100, 1, 1), 2, 3,
                                dimnames = list(c("n", "k"), c("Min", "Max", "Integer")))),

    # Hypotheses matrix
    shinyMatrix::matrixInput("Hypnums", class = "numeric",
                             cols = list(names = TRUE),
                             rows = list(extend = TRUE, names = TRUE, editableNames = TRUE),
                             value =  matrix(c(0.3, 0, 0.05, 0.05, 0.95, 0.95), ncol = 3,
                                             dimnames = list(letters[1:2], c("mu", "var_u", "var_e")))),

    # Constraint matrix
    shinyMatrix::matrixInput("ConMat", class = "numeric",
                             cols = list(names = TRUE),
                             rows = list(names = TRUE),
                             value =  matrix(c(0.2, rep(NA, 5)), ncol = 3,
                                             dimnames = list(letters[1:2], c("s", "n", "k")))),

    # Objectives matrix
    shinyMatrix::matrixInput("ObMat", class = "numeric",
                             cols = list(names = TRUE),
                             rows = list(names = TRUE),
                             value =  matrix(c(NA, NA, 2, NA, 5, NA), ncol = 3,
                                             dimnames = list(letters[1:2], c("s", "n", "k")))),

    shiny::numericInput("test", "test", value = 20),

    # number of initial DoE
    shiny::numericInput("size", "Inital DoE size", value = 20),

    # number of MC evals
    shiny::numericInput("N", "Number of MC evals", value = 100),

    # Button to evaluate
    shiny::actionButton("initButton", "Initialise DoE"),

    shiny::tableOutput("table"),

    # Button to evaluate
    shiny::actionButton("solveButton", "Get approximation set"),

    shiny::tableOutput("tableb"),
  )

  server <- function(input, output, session) {

    shiny::observeEvent(input$Hypnums, {
      shinyMatrix::updateMatrixInput(session, inputId = "ConMat",
                                      value =  matrix(input$ConMat, ncol = 3,
                                                      dimnames = list(rownames(input$Hypnums)[1:2],
                                                                      c("s", "n", "k"))))
      shinyMatrix::updateMatrixInput(session, inputId = "ObMat",
                                     value =  matrix(input$ObMat, ncol = 3,
                                                     dimnames = list(rownames(input$Hypnums)[1:2],
                                                                     c("s", "n", "k"))))
    })

    ds <- shiny::reactive(
      data.frame(name = rownames(input$DSnums),
                              low = input$DSnums[,1],
                              up = input$DSnums[,2],
                              int = input$DSnums[,3])
    )

    hyps <- shiny::reactive({
      m <- data.frame(input$Hypnums)
      m <- m[rowSums(is.na(m)) != ncol(m),]
      m <- m[, colSums(is.na(m)) != nrow(m)]
      m
    })

    get_cons <- shiny::reactive({
      # For each non-NA cell, need a row in the constraint data frame
      m <- input$ConMat
      cons <- NULL
      for(i in 1:nrow(m)){
        for(j in 1:ncol(m)){
          if(!is.na(m[i,j])) cons <- rbind(cons, c(j, i, m[i,j]))
        }
      }

      if(!is.null(cons)){
        cons <- data.frame(name = letters[1:nrow(cons)],
                           out_i = cons[,1],
                           hyp_i = cons[,2],
                           nom = cons[,3],
                           delta = rep(0.975, nrow(cons)),
                           stoch = rep(TRUE, nrow(cons)))
      }

      cons
    })

    get_ob <- shiny::reactive({
      # For each non-NA cell, need new rows in the objectives data frame
      m <- input$ObMat
      ob <- NULL
      for(i in 1:nrow(m)){
        for(j in 1:ncol(m)){
          if(!is.na(m[i,j])) ob <- rbind(ob, c(j, i, m[i,j]))
        }
      }

      if(!is.null(ob)){
        ob <- data.frame(name = paste0("f", 1:nrow(ob)),
                           out_i = ob[,1],
                           hyp_i = ob[,2],
                           weight = ob[,3],
                           stoch = rep(FALSE, nrow(ob)))
        ob$weight <- ob$weight/sum(ob$weight)
        ob$name <- as.character(ob$name)
      }

      ob
    })

    rv <- shiny::reactiveValues(DoE=NULL,count=0)

    #make_DoE <-
    shiny::observeEvent(input$initButton,{
      consss <- get_cons()
      design_space <- ds()

      hypotheses <- hyps()

      DoE <- init_DoE(input$size, design_space)

      N <- 100
      DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=input$N, sim=sim_trial)))
      DoE$N <- input$N
      rv$DoE <- DoE
    })

    b <- shiny::eventReactive(input$solveButton,{
      design_space <- ds()
      #make_DoE()
      DoE <- rv$DoE
      objectives <- get_ob()
      constraints <- get_cons()

      to_model <- unique.data.frame(rbind(objectives[, c("out_i", "hyp_i")],
                                          constraints[, c("out_i", "hyp_i")]))

      models <- fit_models(DoE, to_model, design_space)

      b <- best(design_space, models, DoE, objectives, constraints, to_model)
      b
    })

    output$table <- shiny::renderTable({
      rv$DoE
      #get_cons()
    })

    output$tableb <- shiny::renderTable({
      b()
      #get_ob()
    })
  }

  shiny::shinyApp(ui, server, ...)
}
