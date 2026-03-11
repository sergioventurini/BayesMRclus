# Functions
logpost_beta <- function(beta, gamma, prior, data, log = TRUE) {
  # unpack data
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  # unpack prior
  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  # precompute per-SNP constants
  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  # ensure beta is a numeric vector
  beta <- as.numeric(beta)
  n_beta <- length(beta)
  n_snps <- length(gammahat_j)

  # replicate gammahat_j and Gammahat_j across beta values
  # result: n_snps x n_beta matrices
  gammahat_mat <- matrix(rep(gammahat_j, n_beta), ncol = n_beta)
  Gammahat_mat <- matrix(rep(Gammahat_j, n_beta), ncol = n_beta)
  psi2_j_mat <- matrix(rep(psi2_j, n_beta), ncol = n_beta)
  tau2_j_mat <- matrix(rep(tau2_j, n_beta), ncol = n_beta)

  # beta vector replicated as row vector for broadcasting
  beta_mat <- matrix(rep(beta, each = n_snps), nrow = n_snps)

  a_j <- (beta_mat^2) * psi2 + tau2_j_mat    # n_snps x n_beta
  c_beta <- beta_mat * psi2
  v_j <- a_j * psi2_j_mat - c_beta^2

  # likelihood contributions
  ll_g <- a_j * (gammahat_mat - gamma)^2 / v_j
  ll_G <- psi2_j_mat * (Gammahat_mat - beta_mat * gamma)^2 / v_j
  ll_gG <- -2 * c_beta * (gammahat_mat - gamma) * (Gammahat_mat - beta_mat * gamma) / v_j

  # total loglikelihood per beta
  loglik <- colSums(log(v_j) + ll_g + ll_G + ll_gG)

  # prior contribution
  logprior_beta <- (beta - mu_beta)^2 / sigma2_beta

  # log-posterior
  res <- -0.5 * (loglik + logprior_beta)

  if (!log) res <- exp(res)

  res
}

bayesmr_prior <- function(gammaj = list(psi2 = 1), Gammaj = list(tau2 = 1),
                          gamma = list(mean = 0, var = 1),
                          beta = list(mean = 0, var = 1)){
  prior <- list()
  for (arg in names(formals(sys.function())))
    prior[[arg]] <- get(arg)
  prior
}

gamma_beta_post <- function(gamma, beta, data, prior, log = TRUE) {
  # unpack data
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  # unpack prior
  mu_gamma <- prior[["gamma"]][["mean"]]
  sigma2_gamma <- prior[["gamma"]][["var"]]
  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  # precompute per-SNP constants
  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  gamma_len <- length(gamma)
  beta_len <- length(beta)
  res <- matrix(NA, nrow = gamma_len, ncol = beta_len)

  for (g in 1:gamma_len) {
    p_gamma <- (gamma[g] - mu_gamma)^2/sigma2_gamma
    for (b in 1:beta_len) {
      p_beta <- (beta[b] - mu_beta)^2/sigma2_beta
      a_j <- beta[b]^2*psi2 + tau2_j
      c_beta <- beta[b]*psi2
      v_j <- psi2_j*a_j - c_beta^2  # same as beta[b]^2*psi2*sigma2_X + psi2_j*tau2_j
      p_gj <- a_j*(gammahat_j - gamma[g])^2
      p_Gj <- psi2_j*(Gammahat_j - beta[b]*gamma[g])^2
      p_gj_Gj <- -2*c_beta*(gammahat_j - gamma[g])*(Gammahat_j - beta[b]*gamma[g])
      
      res[g, b] <- -0.5*(sum(log(v_j) + (p_gj + p_Gj + p_gj_Gj)/v_j) + p_gamma + p_beta)
    }
  }
  res <- res[, , drop = TRUE]

  if (!log) res <- exp(res)

  res
}

gamma_beta_logpost_optim <- function(par, data, prior) {
  gamma <- par[1]
  beta  <- par[2]

  as.numeric(
    gamma_beta_post(
      gamma = gamma,
      beta  = beta,
      data  = data,
      prior = prior,
      log   = TRUE
    )
  )
}

bayesmr_noclus_optim <- function(
  data,
  prior,
  init = c(gamma = 0, beta = 0),
  method = "BFGS",
  control = list(fnscale = -1, reltol = 1e-10, maxit = 1000)) {
  opt <- optim(
    par = init,
    fn = gamma_beta_logpost_optim,
    data = data,
    prior = prior,
    method = method,
    hessian = TRUE,
    control = control
  )

  vcov <- tryCatch(
    solve(-opt$hessian),
    error = function(e) matrix(NA, 2, 2)
  )

  # return results
  list(
    par = setNames(opt$par, c("gamma", "beta")),
    value = opt$value,
    convergence = opt$convergence,
    hessian = opt$hessian,
    vcov = vcov,
    sd = sqrt(diag(vcov)),
    optim = opt
  )
}

# Install packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("shiny", quietly = TRUE)) {
  install.packages("shiny")
}
if (!requireNamespace("mr.raps", quietly = TRUE)) {
  devtools::install_github("qingyuanzhao/mr.raps")
}
library(ggplot2)
library(shiny)
library(mr.raps)

# Slider with a specified column width
slider_cols <- function(cols, inputId, label, min, max, value, step = NULL) {
  column(cols,
    tags$div(withMathJax(label), style = "text-align: left; margin-bottom: 5px;"),
    sliderInput(inputId, NULL, min, max, value, step = step)
  )
}

# UI
ui <- fluidPage(
  titlePanel("Interactive (Unnormalized) Joint Posterior"),

  withMathJax(),

  tags$style(HTML("
    .form-group {
      text-align: left !important;
    }
    .form-group > label {
      display: block;
      text-align: left !important;
    }
  ")),

  fluidRow(
    column(3, 
      selectInput(
        inputId = "dataset",          # name used in server
        label = "Choose a dataset:",  # label shown to user
        choices = c("bmi.sbp", "bmi.ais", "bmi.cad", "crp.cad", "simulated"),  # available options
        selected = "bmi.sbp"             # default selection
      )
    ),
    column(3,
      tags$b("P-value for SNPs selection:"),
      textInput("pval", NULL, value = "1e-4")
    ),
    column(3,
      tags$b("Number of selected SNPs:"),
      div(
        textOutput("n_selected"),
        class = "form-control", style = "background-color: #f9f9f9;"
      )
    ),
    column(3,
      conditionalPanel(
        condition = "input.dataset == 'simulated'",
        numericInput("sim_seed", "Random seed:", value = 123, step = 1)
      )
    )
  ),

  fluidRow(column(12, tags$h4("Axis Ranges", style = "margin-top: 20px;"))),
  fluidRow(
    column(3, tags$div("$$\\gamma \\text{ min}$$"), textInput("gamma_min", NULL, value = "-0.05")),
    column(3, tags$div("$$\\gamma \\text{ max}$$"), textInput("gamma_max", NULL, value = "0.05")),
    column(3, tags$div("$$\\beta \\text{ min}$$"), textInput("beta_min", NULL, value = "-1")),
    column(3, tags$div("$$\\beta \\text{ max}$$"), textInput("beta_max", NULL, value = "1"))
  ),

  fluidRow(column(12, tags$h4("Prior Parameters", style = "margin-top: 30px;"))),

  fluidRow(
    slider_cols(4, "mu_gamma", "$$\\mu_\\gamma$$", -2, 2, 0, step = 0.01),
    slider_cols(4, "mu_beta", "$$\\mu_\\beta$$", -2, 2, 0, step = 0.01),
    slider_cols(4, "psi2", "$$\\psi^2$$", 0, 10, 0.0001, step = 0.0001)
  ),
  fluidRow(
    slider_cols(4, "sigma2_gamma", "$$\\sigma^2_\\gamma$$", 0, 100, 100, step = 0.01),
    slider_cols(4, "sigma2_beta", "$$\\sigma^2_\\beta$$", 0, 100, 100, step = 0.01),
    slider_cols(4, "tau2", "$$\\tau^2$$", 0, 10, 0.0001, step = 0.0001)
  ),

  fluidRow(
    column(4, tags$h4("Plot Settings", style = "margin-top: 30px;")),
    column(4, selectInput("resolution", "Grid Resolution",
      choices = c("Low" = 50, "Medium" = 100, "High" = 150), selected = 50)),
    column(4,
        tags$div(
          style = "border: 1px solid #ddd; padding: 10px; background-color: #f9f9f9; border-radius: 5px;",
          "Note: the coloured dot and solid lines represent the joint posterior global maximum."))
  ),

  fluidRow(column(12, plotOutput("joint_plot", height = "700px"))),

  fluidRow(
    slider_cols(6, "beta_for_gamma", "$$\\beta \\text{ value}$$", -3, 3, 0, step = 0.01),
    slider_cols(6, "gamma_for_beta", "$$\\gamma \\text{ value}$$", -3, 3, 0, step = 0.01)
  ),

  fluidRow(
    column(6, plotOutput("gamma_plot", height = "300px")),
    column(6, plotOutput("beta_plot", height = "300px"))
  )
)

# Server
server <- function(input, output, session) {

  dataset_map <- list(
    bmi.sbp       = "mr.raps",
    bmi.ais       = "mr.raps",
    bmi.cad       = "mr.raps",
    crp.cad       = "mr.raps"
  )

  datasetInput <- reactive({
    if (input$dataset == "simulated") {
      set.seed(input$sim_seed)
      nsnips <- 1500
      g <- 0.2
      b <- 0.4
      psi2 <- 0.1
      tau2 <- 0.1
      sigma2_X <- rep(0.0001, nsnips)
      sigma2_Y <- rep(0.00001, nsnips)
      gamma_j <- rnorm(nsnips, mean = g, sd = sqrt(psi2))
      Gamma_j <- rnorm(nsnips, mean = b*gamma_j, sd = sqrt(tau2))
      gammahat_j <- rnorm(nsnips, mean = gamma_j, sd = sqrt(sigma2_X))
      Gammahat_j <- rnorm(nsnips, mean = Gamma_j, sd = sqrt(sigma2_Y))
      data.frame(
        beta.exposure = gammahat_j,
        beta.outcome  = Gammahat_j,
        se.exposure   = sqrt(sigma2_X),
        se.outcome    = sqrt(sigma2_Y),
        pval.selection = rbeta(nsnips, shape1 = 1e-1, shape2 = 1e4)
      )
    } else {
      pkg <- dataset_map[[input$dataset]]
      get(input$dataset, paste0("package:", pkg))
    }
  })

  fixed_var <- "pval.selection"

  filteredData <- reactive({
    data <- datasetInput()
    if (fixed_var %in% names(data)) {
      pval <- suppressWarnings(as.numeric(input$pval))
      if (is.na(pval) || pval <= 0 || pval >= 1) {
        return(data)
      } else {
        return(subset(data, data[[fixed_var]] <= pval))
      }
    } else {
      data
    }
  })
  
  output$n_selected <- renderText({
    nrow(filteredData())
  })

  plotParams <- reactive({
    list(
      mu_gamma     = input$mu_gamma,
      mu_beta      = input$mu_beta,
      sigma2_gamma = input$sigma2_gamma,
      sigma2_beta  = input$sigma2_beta,
      psi2         = input$psi2,
      tau2         = input$tau2,
      gamma_val    = input$gamma_for_beta,
      beta_val     = input$beta_for_gamma
    )
  })

  output$joint_plot <- renderPlot({
    gamma_min <- suppressWarnings(as.numeric(input$gamma_min))
    gamma_max <- suppressWarnings(as.numeric(input$gamma_max))
    beta_min  <- suppressWarnings(as.numeric(input$beta_min))
    beta_max  <- suppressWarnings(as.numeric(input$beta_max))

    if (is.na(gamma_min) || is.na(gamma_max) || gamma_min >= gamma_max) {
      gamma_min <- -0.05; gamma_max <- 0.05
    }
    if (is.na(beta_min) || is.na(beta_max) || beta_min >= beta_max) {
      beta_min <- -1; beta_max <- 1
    }

    data <- data.frame(
      beta_exposure = filteredData()[, "beta.exposure"],
      beta_outcome = filteredData()[, "beta.outcome"],
      se_exposure = filteredData()[, "se.exposure"],
      se_outcome = filteredData()[, "se.outcome"]
    )

    p <- plotParams()

    res <- as.numeric(input$resolution)
    gamma_vals <- seq(gamma_min, gamma_max, length.out = res)
    beta_vals  <- seq(beta_min, beta_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = p$mu_gamma, var = p$sigma2_gamma),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    post_vals <- gamma_beta_post(gamma_vals, beta_vals, data, prior, log = TRUE)

    df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
    df_plot$posterior <- as.vector(post_vals)

    map <- bayesmr_noclus_optim(data, prior, init = rnorm(2, mean = 0, sd = 20))

    ggplot(df_plot, aes(x = gamma, y = beta, z = posterior)) +
      geom_contour_filled(bins = 20) +
      scale_fill_viridis_d(option = "C") +
      geom_vline(xintercept = map$par["gamma"], linetype = "solid", color = "#21908C") +
      geom_hline(yintercept = map$par["beta"], linetype = "solid", color = "#21908C") +
      labs(x = expression(gamma), y = expression(beta), fill = "Posterior") +
      coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
      geom_point(aes(x = map$par["gamma"], y = map$par["beta"]), color = "#21908C", size = 3) +
      theme_minimal(base_size = 14)
  })

  output$gamma_plot <- renderPlot({
    gamma_min <- suppressWarnings(as.numeric(input$gamma_min))
    gamma_max <- suppressWarnings(as.numeric(input$gamma_max))

    if (is.na(gamma_min) || is.na(gamma_max) || gamma_min >= gamma_max) {
      gamma_min <- -0.05; gamma_max <- 0.05
    }

    data <- data.frame(
      beta_exposure = filteredData()[, "beta.exposure"],
      beta_outcome = filteredData()[, "beta.outcome"],
      se_exposure = filteredData()[, "se.exposure"],
      se_outcome = filteredData()[, "se.outcome"]
    )

    p <- plotParams()

    res <- as.numeric(input$resolution)
    gamma_vals <- seq(gamma_min, gamma_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = p$mu_gamma, var = p$sigma2_gamma),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    post_vals <- gamma_beta_post(gamma_vals, p$beta_val, data, prior, log = TRUE)
    df <- data.frame(gamma = gamma_vals, posterior = as.numeric(post_vals))

    map <- bayesmr_noclus_optim(data, prior, init = rnorm(2, mean = 0, sd = 20))

    ggplot(df, aes(x = gamma, y = posterior)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      # geom_vline(xintercept = p$mu_gamma, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = map$par["gamma"], linetype = "dashed", color = "black") +
      labs(x = expression(gamma), y = paste("Marginal Posterior at β =", p$beta_val)) +
      theme_minimal(base_size = 14)
  })

  output$beta_plot <- renderPlot({
    beta_min <- suppressWarnings(as.numeric(input$beta_min))
    beta_max <- suppressWarnings(as.numeric(input$beta_max))

    if (is.na(beta_min) || is.na(beta_max) || beta_min >= beta_max) {
      beta_min <- -1; beta_max <- 1
    }

    data <- data.frame(
      beta_exposure = filteredData()[, "beta.exposure"],
      beta_outcome = filteredData()[, "beta.outcome"],
      se_exposure = filteredData()[, "se.exposure"],
      se_outcome = filteredData()[, "se.outcome"]
    )

    p <- plotParams()

    res <- as.numeric(input$resolution)
    beta_vals <- seq(beta_min, beta_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = p$mu_gamma, var = p$sigma2_gamma),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    # post_vals <- gamma_beta_post(p$gamma_val, beta_vals, data, prior, log = FALSE)
    post_vals <- logpost_beta(beta_vals, p$gamma_val, prior, data, log = TRUE)
    df <- data.frame(beta = beta_vals, posterior = as.numeric(post_vals))

    map <- bayesmr_noclus_optim(data, prior, init = rnorm(2, mean = 0, sd = 20))

    ggplot(df, aes(x = beta, y = posterior)) +
      geom_line(color = "darkorange", linewidth = 1.2) +
      # geom_vline(xintercept = p$mu_beta, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = map$par["beta"], linetype = "dashed", color = "black") +
      labs(x = expression(beta), y = paste("Marginal Posterior at γ =", p$gamma_val)) +
      theme_minimal(base_size = 14)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
