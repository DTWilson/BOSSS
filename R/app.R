BOSSSapp <- function(...) {

  set.seed(1)

  get_det_obj <- function(design)
  {
    o <- matrix(design, ncol = 2)[,1:2]
    c(s = NA, p = o[1], c = o[2])
  }

  det_obj1 <- function(design)
  {
    if(is.null(dim(design))) design <- matrix(design, nrow = 1)
    design[,1]
  }

  det_obj2 <- function(design)
  {
    if(is.null(dim(design))) design <- matrix(design, nrow = 1)
    design[,2]
  }

  ## use GPareto fastfun for deterministic objectives

  ui <- shiny::fluidPage(
    theme = bslib::bs_theme(bootswatch = "united"),

    shiny::titlePanel(
      shiny::h1("BOSSS",
         shiny::h5("Bayesian Optimisation for Sample Size with Simulation"))
    ),

    shiny::tabsetPanel(
      shiny::tabPanel("Problem",
        # Simulation code
        shiny::textAreaInput("simraw", "Simulation code",
                             value = "n <- design[1]; k <- design[2]
        mu <- hypothesis[,1]; var_u <- hypothesis[,2]; var_e <- hypothesis[,3]

        m <- n/k
        s_c <- sqrt(var_u + var_e/m)
        x0 <- stats::rnorm(k, 0, s_c); x1 <- stats::rnorm(k, mu, s_c)
        return(c(s = stats::t.test(x0, x1)$p.value >= 0.05, p = n, c = k))",
                             width = '800px', height = '400px'),

        shiny::textInput("desnames", "Design variable names:",
                         "n,k"),

        shiny::textInput("hypnames", "Model parameter names:",
                         "mu,var_u,var_e"),

        shiny::textInput("outnames", "Output variable names:",
                         "s,p,c"),

        shiny::actionButton("checkSim", "Read in simulation")
      ),

      shiny::tabPanel("Initialisation",
        # Design space matrix
        shinyMatrix::matrixInput("DSnums", label = "Design space", class = "numeric",
                                  cols = list(names = TRUE),
                                  rows = list(names = TRUE),
                                  value =  matrix(c(rep(1,6), ncol=3))),

        # Hypotheses matrix
        shinyMatrix::matrixInput("Hypnums", label = "Hypotheses", class = "numeric",
                                  cols = list(names = TRUE),
                                  rows = list(extend = TRUE, names = TRUE, editableNames = TRUE),
                                  value = matrix(c(rep(1,3), ncol=3))),

        # Constraint matrix
        shinyMatrix::matrixInput("ConMat", label = "Constraints", class = "numeric",
                                 cols = list(names = TRUE),
                                 rows = list(names = TRUE),
                                 value =  matrix(c(rep(1,3), ncol=3))),

        # Objectives matrix
        shinyMatrix::matrixInput("ObMat", label = "Objectives", class = "numeric",
                                 cols = list(names = TRUE),
                                 rows = list(names = TRUE),
                                 value =  matrix(c(rep(1,3), ncol=3))),

        # number of initial DoE
        shiny::numericInput("size", "Inital DoE size", value = 20),

        # number of MC evals
        shiny::numericInput("N", "Number of MC evals", value = 100),

        # Button to evaluate
        shiny::actionButton("initButton", "Initialise DoE"),

        shiny::tableOutput("table")
      ),

      shiny::tabPanel("Iteration",
        # Button to iterate
        shiny::actionButton("iterButton", "Perform an iteration"),

        shiny::tableOutput("tableb")#,

        #shiny::plotOutput("ASgraph"),

        #shiny::plotOutput("Trajgraph")
      )
    )
  )

  server <- function(input, output, session) {
    thematic::thematic_shiny()

    shiny::observeEvent(input$Hypnums, {

      # Set up constraint and objective matrices
      out_name_split <- strsplit(input$outnames, ",")[[1]]
      hyp_num <- nrow(input$Hypnums)

      shinyMatrix::updateMatrixInput(session, inputId = "ConMat",
                                      value =  matrix(rep(NA, length(out_name_split)*hyp_num),
                                                      ncol = length(out_name_split),
                                                      dimnames = list(rownames(input$Hypnums),
                                                                      out_name_split)))
      shinyMatrix::updateMatrixInput(session, inputId = "ObMat",
                                     value =  matrix(rep(NA, length(out_name_split)*hyp_num),
                                                     ncol = length(out_name_split),
                                                     dimnames = list(rownames(input$Hypnums),
                                                                     out_name_split)))
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
                           stoch = rep(TRUE, nrow(ob)))
        ob$weight <- ob$weight/sum(ob$weight)
        ob$name <- as.character(ob$name)
      }

      ob
    })

    rv <- shiny::reactiveValues(DoE=NULL, PS=NULL, models=NULL, traj=NULL, sim_trial = NULL)

    shiny::observeEvent(input$checkSim,{
      # Read in the simulation and turn into a function
      rv$sim_trial <- eval(parse(text = paste('f <- function(design, hypothesis) {', input$simraw, '}', sep='')))

      # Set up design space matrix
      ds_name_split <- strsplit(input$desnames, ",")[[1]]

      shinyMatrix::updateMatrixInput(session, inputId = "DSnums",
                                     value =  matrix(c(10,10,500,50,1,1),#rep(NA, 3*length(ds_name_split)),
                                                     ncol = 3,
                                                     dimnames = list(ds_name_split,
                                                                     c("Min", "Max", "Integer"))))

      # Set up hypothesis matrix
      hyp_name_split <- strsplit(input$hypnames, ",")[[1]]

      shinyMatrix::updateMatrixInput(session, inputId = "Hypnums",
                                     value =  matrix(c(0,0.05,0.95),#rep(NA, length(hyp_name_split)),
                                                     ncol = length(hyp_name_split),
                                                     dimnames = list("a", hyp_name_split)))
    })

    shiny::observeEvent(input$initButton,{

      shiny::withProgress(message = 'Initialising', value = 0, {

        design_space <- ds()
        objectives <- get_ob()
        constraints <- get_cons()
        hypotheses <- hyps()

        out_dim <- ncol(input$ObMat)

        DoE <- init_DoE(input$size, design_space)

        shiny::incProgress(1/3, detail = "Running simulations")

        DoE <- cbind(DoE, t(apply(DoE, 1, calc_rates, hypotheses=hypotheses, N=input$N, sim=rv$sim_trial)))
        DoE$N <- input$N
        rv$DoE <- DoE

        to_model <- rbind(constraints[,c("out_i", "hyp_i")], objectives[,c("out_i", "hyp_i")])
        to_model <- unique(to_model)

        incProgress(2/3, detail = "Fitting models")

        rv$models <- fit_models(DoE, to_model, design_space, objectives, out_dim)

        pf_out <- pareto_front(design_space, rv$models, rv$DoE, objectives, constraints, to_model, out_dim, get_det_obj)
        rv$PS <- pf_out[[1]][,1:nrow(objectives)]

      })
    })

    shiny::observeEvent(input$iterButton,{
      # choose a random point in the design space
      design_space <- ds()
      hypotheses <- hyps()
      objectives <- get_ob()
      constraints <- get_cons()

      out_dim <- ncol(input$ObMat)

      to_model <- rbind(constraints[,c("out_i", "hyp_i")], objectives[,c("out_i", "hyp_i")])
      to_model <- unique(to_model)

      shiny::withProgress(message = 'Iterating', value = 0, {

        incProgress(1/4, detail = "Optimising")

        opt <- RcppDE::DEoptim(ehi_infill, lower=design_space$low, upper=design_space$up,
                               control=list(trace=FALSE, itermax=100, reltol=1e-1, steptol=50),
                               N=input$N, pf=rv$PS, mod=rv$models, design_space=design_space, constraints=constraints,
                               objectives=objectives, det_obj=get_det_obj, out_dim=out_dim, to_model=to_model)

        sol <- as.numeric(opt$optim$bestmem)

        incProgress(2/4, detail = "Evaluating selected design")

        #sol[1:2] <- round(sol[1:2])
        y <- calc_rates(sol, hypotheses=hypotheses, N=input$N, sim=rv$sim_trial)

        rv$DoE <- rbind(rv$DoE, c(sol, y, input$N))

        incProgress(3/4, detail = "Updating models")

        rv$models <- fit_models(rv$DoE, to_model, design_space, objectives, out_dim)

        pf_out <- pareto_front(design_space, rv$models, rv$DoE, objectives, constraints, to_model, out_dim, get_det_obj)
        rv$PS <- pf_out[[1]][,1:nrow(objectives)]

        ref <- design_space$up*objectives$weight

        current <- emoa::dominated_hypervolume(t(rv$PS), ref)
        rv$traj <- rbind(rv$traj, c(-opt$optim$bestval, current))

      })
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
        ggplot2::geom_point() #data=df[2:(nrow(df)-1),])
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
