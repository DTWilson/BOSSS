BOSSSapp <- function(...) {

  set.seed(1)

  sim_trial <- function(design, hypothesis)
  {
    n <- design[1]; k <- design[2]
    mu <- hypothesis[,1]; var_u <- hypothesis[,2]; var_e <- hypothesis[,3]

    m <- n/k
    s_c <- sqrt(var_u + var_e/m)
    x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
    c(s = stats::t.test(x0, x1)$p.value >= 0.05, p = n, c = k)
  }

  get_det_obj <- function(design)
  {
    o <- matrix(design, ncol = 2)[,1:2]
    c(s = NA, p = o[1], c = o[2])
  }

  #num_hyp <- 2

  ui <- shiny::fluidPage(

    # Simulation code
    shiny::textAreaInput("simraw", "Simulation code",
                         value = "n <- design[1]; k <- design[2]
    mu <- hypothesis[,1]; var_u <- hypothesis[,2]; var_e <- hypothesis[,3]

    m <- n/k
    s_c <- sqrt(var_u + var_e/m)
    x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
    c(s = stats::t.test(x0, x1)$p.value >= 0.05, p = n, c = k)"),

    # Deterministic objective code
    shiny::textAreaInput("objraw", "Deterministic objective code"),

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

    # number of initial DoE
    shiny::numericInput("size", "Inital DoE size", value = 20),

    # number of MC evals
    shiny::numericInput("N", "Number of MC evals", value = 100),

    shiny::actionButton("checkSim", "Read in simulation"),

    # Button to evaluate
    shiny::actionButton("initButton", "Initialise DoE"),

    shiny::tableOutput("table"),

    shiny::tableOutput("tableb"),

    shiny::plotOutput("ASgraph"),

    shiny::plotOutput("Trajgraph"),

    # Button to iterate
    shiny::actionButton("iterButton", "Perform an iteration"),
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

    rv <- shiny::reactiveValues(DoE=NULL, PS=NULL, models=NULL, traj=NULL, sim_trial = NULL)

    shiny::observeEvent(input$checkSim,{
      rv$sim_trial <- eval(parse(text = paste('f <- function(design, hypothesis) {', input$simraw, '}', sep='')))
    })

    shiny::observeEvent(input$initButton,{
      design_space <- ds()
      objectives <- get_ob()
      constraints <- get_cons()
      hypotheses <- hyps()

      DoE <- init_DoE(input$size, design_space)

      DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=input$N, sim=rv$sim_trial)))
      DoE$N <- input$N
      rv$DoE <- DoE

      to_model <- data.frame(out_i = c(1),
                             hyp_i = c(1))

      rv$models <- fit_models(DoE, to_model, design_space)

      rv$PS <- best(design_space, rv$models, rv$DoE, objectives, constraints, to_model, get_det_obj)
    })

    shiny::observeEvent(input$iterButton,{
      # choose a random point in the design space
      design_space <- ds()
      hypotheses <- hyps()
      objectives <- get_ob()
      constraints <- get_cons()

      to_model <- data.frame(out_i = c(1),
                             hyp_i = c(1))

      opt <- RcppDE::DEoptim(exp_improve, lower=design_space$low, upper=design_space$up,
                             control=list(trace=FALSE),
                             N=input$N, PS=rv$PS, mod=rv$models, design_space=design_space, constraints=constraints,
                             objectives=objectives, get_det_obj=get_det_obj, out_dim=3, to_model = to_model)

      sol <- as.numeric(opt$optim$bestmem)

      sol[1:2] <- round(sol[1:2])
      y <- calc_rates(sol, hypotheses=hypotheses, N=input$N, sim=rv$sim_trial)

      rv$DoE <- rbind(rv$DoE, c(sol, y, input$N))

      #to_model <- unique.data.frame(rbind(objectives[, c("out_i", "hyp_i")],
      #                                    constraints[, c("out_i", "hyp_i")]))
      to_model <- data.frame(out_i = c(1),
                             hyp_i = c(1))

      rv$models <- fit_models(rv$DoE, to_model, design_space)

      rv$PS <- best(design_space, rv$models, rv$DoE, objectives, constraints, to_model, get_det_obj)

      ref <- design_space$up*objectives$weight
      PS2 <- as.matrix(rv$PS[, objectives$name])
      current <- mco::dominatedHypervolume(PS2, ref)

      rv$traj <- rbind(rv$traj, c(-opt$optim$bestval, current))
    })

    output$table <- shiny::renderTable({
      rv$DoE
    })

    output$tableb <- shiny::renderTable({
      rv$PS
    })

    output$ASgraph <- shiny::renderPlot({
      objectives <- get_ob()
      df <- as.data.frame(rv$PS[,objectives$name])
      ## Extend to include extreme points for plotting
      #df <- rbind(c(min(df[,1]), 100*objectives$weight[2]),
      #            df,
      #            c(500*objectives$weight[1], min(df[,2])))

      ggplot2::ggplot(df, ggplot2::aes(f1, f2)) +
        ggplot2::geom_step(linetype = 2) +
        ggplot2::geom_point(data=df[2:(nrow(df)-1),])
    })

    output$Trajgraph <- shiny::renderPlot({
      traj <- rv$traj
      if(!is.null(traj)){

        ## scale EI trajectory to have same min and max as DH
        traj[,1] <- (max(traj[,2]) - min(traj[,2]))*(traj[,1] - min(traj[,1]))/(max(traj[,1] - min(traj[,1]))) + min(traj[,2])

        n <- nrow(traj)
        df <- data.frame(v = c(traj[,2], traj[,1]),
                         t = rep(c("DH", "EI"), each = n),
                         i = rep(1:n, 2))
        ggplot2::ggplot(df, ggplot2::aes(i, v, colour = t)) +
          ggplot2::geom_point() +
          ggplot2::geom_line()
      }
    })
  }

  shiny::shinyApp(ui, server, ...)
}
