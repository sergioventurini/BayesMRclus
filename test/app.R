# Functions
bayesmr_prior <- function(gammaj = list(psi2 = 1), Gammaj = list(tau2 = 1),
                          gamma = list(mean = 0, var = 1),
                          beta = list(mean = 0, var = 1)){
  prior <- list()
  for (arg in names(formals(sys.function())))
    prior[[arg]] <- get(arg)
  prior
}

gamma_beta_post <- function(gamma, beta, data, prior, log = TRUE, verbose = TRUE) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  mu_gamma <- prior[["gamma"]][["mean"]]
  sigma2_gamma <- prior[["gamma"]][["var"]]
  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  gamma_len <- length(gamma)
  beta_len <- length(beta)
  res <- matrix(NA, nrow = gamma_len, ncol = beta_len)

  for (g in 1:gamma_len) {
    p_gamma <- (gamma[g] - mu_gamma)^2/sigma2_gamma
    p_gj <- (gammahat_j - gamma[g])^2/psi2_j
    for (b in 1:beta_len) {
      if (verbose)
        message("Computing gamma/beta joint posterior for gamma = ", g, " and beta = ", b)

      p_beta <- (beta[b] - mu_beta)^2/sigma2_beta
      p_Gj <- (Gammahat_j - beta[b]*gamma[g])^2/(beta[b]^2*psi2 + tau2_j)
      
      res[g, b] <- -0.5*(sum(p_gj) + sum(p_Gj) + p_gamma + p_beta)
    }
  }

  if (!log)
    res <- exp(res)

  res
}

beta_marg_post_drv <- function(beta, gamma, data, prior, beta_min = -20, beta_max = 20) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  tau2_j <- sigma2_Y + tau2

  res <- numeric(length(beta))
  for (b in 1:length(beta)) {
    h_j_beta <- beta[b]^2*psi2 + tau2_j

    res[b] <- beta[b]^2*(gamma*psi2*sum(Gammahat_j/h_j_beta^2)) +
              beta[b]*(gamma^2*sum(tau2_j/h_j_beta^2) - psi2*sum(Gammahat_j^2/h_j_beta^2) + 1/sigma2_beta) -
              (gamma*sum(Gammahat_j*tau2_j/h_j_beta^2) + mu_beta/sigma2_beta)
  }

  res
}

gamma_post_mode <- function(beta, data, prior) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  mu_gamma <- prior[["gamma"]][["mean"]]
  sigma2_gamma <- prior[["gamma"]][["var"]]
  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2
  h2_j <- (beta^2)*psi2 + tau2_j
  A_beta <- sum(1/psi2_j) + (beta^2)*sum(1/h2_j) + 1/sigma2_gamma
  B_beta <- sum(gammahat_j/psi2_j) + beta*sum(Gammahat_j/h2_j) + mu_gamma/sigma2_gamma
  res <- B_beta/A_beta

  res
}

find_all_roots <- function(f, ..., lower, upper, n = 1000,
                           tol_x = 1e-10, tol_f = 1e-12,
                           eps_small = 1e-8) {
  # Create coarse sampling grid
  xs <- seq(lower, upper, length.out = n)
  vals <- sapply(xs, f, ...)
  
  brackets <- list()
  near_zero_intervals <- list()
  
  for (i in seq_len(n - 1)) {
    # Exact zero at sample point
    if (vals[i] == 0) {
      brackets[[length(brackets) + 1]] <- c(xs[i], xs[i])
    }
    # Sign change
    else if (vals[i] * vals[i + 1] < 0) {
      brackets[[length(brackets) + 1]] <- c(xs[i], xs[i + 1])
    }
    # Possible multiple root (no sign change but small values)
    else if (min(abs(vals[i]), abs(vals[i + 1])) < eps_small) {
      near_zero_intervals[[length(near_zero_intervals) + 1]] <- c(xs[i], xs[i + 1])
    }
  }
  
  # Function to refine root in bracket with uniroot
  refine_root <- function(a, b) {
    if (a == b) return(a)  # already exact
    tryCatch({
      r <- uniroot(f, ..., lower = a, upper = b, tol = tol_x)$root
      return(r)
    }, error = function(e) NA)
  }
  
  roots <- numeric(0)
  
  # Handle sign-change brackets
  for (br in brackets) {
    roots <- c(roots, refine_root(br[1], br[2]))
  }
  
  # Handle near-zero intervals (possible multiple roots)
  for (br in near_zero_intervals) {
    mid <- mean(br)
    if (abs(f(mid, ...)) < tol_f * 100) {
      # Try small bracket around midpoint
      r <- refine_root(br[1], br[2])
      roots <- c(roots, r)
    }
  }
  
  # Remove NA and duplicates
  roots <- roots[!is.na(roots)]
  roots <- sort(unique(roots))
  
  roots
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

# Prepare data
data("bmi.sbp", package = "mr.raps")
data <- data.frame(
  beta_exposure = bmi.sbp[, "beta.exposure"],
  beta_outcome = bmi.sbp[, "beta.outcome"],
  se_exposure = bmi.sbp[, "se.exposure"],
  se_outcome = bmi.sbp[, "se.outcome"]
)

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

  fluidRow(column(12, tags$h4("Axis Ranges", style = "margin-top: 20px;"))),
  fluidRow(
    column(3, tags$div("$$\\gamma \\text{ min}$$"), textInput("gamma_min", NULL, value = "-0.05")),
    column(3, tags$div("$$\\gamma \\text{ max}$$"), textInput("gamma_max", NULL, value = "0.05")),
    column(3, tags$div("$$\\beta \\text{ min}$$"), textInput("beta_min", NULL, value = "-1")),
    column(3, tags$div("$$\\beta \\text{ max}$$"), textInput("beta_max", NULL, value = "1"))
  ),

  fluidRow(column(12, tags$h4("Prior Parameters", style = "margin-top: 30px;"))),

  fluidRow(
    slider_cols(4, "mu_gamma", "$$\\mu_\\gamma$$", -5, 5, 0, step = 0.1),
    slider_cols(4, "mu_beta", "$$\\mu_\\beta$$", -5, 5, 0, step = 0.1),
    slider_cols(4, "psi2", "$$\\psi^2$$", 0, 10, 0.1, step = 0.01)
  ),
  fluidRow(
    slider_cols(4, "sigma2_gamma", "$$\\sigma^2_\\gamma$$", 0, 10, 0.1, step = 0.01),
    slider_cols(4, "sigma2_beta", "$$\\sigma^2_\\beta$$", 0, 10, 0.1, step = 0.01),
    slider_cols(4, "tau2", "$$\\tau^2$$", 0, 10, 0.1, step = 0.01)
  ),

  fluidRow(
    column(4, tags$h4("Plot Settings", style = "margin-top: 30px;")),
    column(4, selectInput("resolution", "Grid Resolution",
      choices = c("Low" = 50, "Medium" = 100, "High" = 150), selected = 50))
  ),

  fluidRow(column(12, plotOutput("joint_plot", height = "700px"))),

  fluidRow(
    slider_cols(6, "beta_for_gamma", "$$\\beta \\text{ value}$$", -20, 20, 0, step = 0.01),
    slider_cols(6, "gamma_for_beta", "$$\\gamma \\text{ value}$$", -20, 20, 0, step = 0.01)
  ),

  fluidRow(
    column(6, plotOutput("gamma_plot", height = "300px")),
    column(6, plotOutput("beta_plot", height = "300px"))
  )
)

# Server
server <- function(input, output, session) {

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

    post_vals <- gamma_beta_post(gamma_vals, beta_vals, data, prior, log = FALSE, verbose = FALSE)

    df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
    df_plot$posterior <- as.vector(post_vals)

    # find beta modes
    beta_modes <- find_all_roots(beta_marg_post_drv,
                                 lower = beta_min, upper = beta_max,
                                 gamma = p$gamma_val, data = data, prior = prior,
                                 n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                                 eps_small = 1e-8)
    if (length(beta_modes) > 1) beta_modes <- beta_modes[-2]

    # find gamma mode
    gamma_mode <- gamma_post_mode(p$beta_val, data, prior)

    ggplot(df_plot, aes(x = gamma, y = beta, z = posterior)) +
      geom_contour_filled(bins = 20) +
      scale_fill_viridis_d(option = "C") +
      geom_vline(xintercept = p$mu_gamma, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = gamma_mode, linetype = "dashed", color = "black") +
      geom_hline(yintercept = p$mu_beta, linetype = "dashed", color = "gray") +
      geom_hline(yintercept = beta_modes, linetype = "dashed", color = "black") +
      labs(x = expression(gamma), y = expression(beta), fill = "Posterior") +
      coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
      theme_minimal(base_size = 14)
  })

  output$gamma_plot <- renderPlot({
    gamma_min <- suppressWarnings(as.numeric(input$gamma_min))
    gamma_max <- suppressWarnings(as.numeric(input$gamma_max))

    if (is.na(gamma_min) || is.na(gamma_max) || gamma_min >= gamma_max) {
      gamma_min <- -0.05; gamma_max <- 0.05
    }

    p <- plotParams()

    res <- as.numeric(input$resolution)
    gamma_vals <- seq(gamma_min, gamma_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = p$mu_gamma, var = p$sigma2_gamma),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    post_vals <- gamma_beta_post(gamma_vals, p$beta_val, data, prior, log = FALSE, verbose = FALSE)
    df <- data.frame(gamma = gamma_vals, posterior = as.numeric(post_vals))

    # find gamma mode
    gamma_mode <- gamma_post_mode(p$beta_val, data, prior)

    ggplot(df, aes(x = gamma, y = posterior)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_vline(xintercept = p$mu_gamma, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = gamma_mode, linetype = "dashed", color = "black") +
      labs(x = expression(gamma), y = paste("Marginal Posterior at β =", p$beta_val)) +
      theme_minimal(base_size = 14)
  })

  output$beta_plot <- renderPlot({
    beta_min <- suppressWarnings(as.numeric(input$beta_min))
    beta_max <- suppressWarnings(as.numeric(input$beta_max))

    if (is.na(beta_min) || is.na(beta_max) || beta_min >= beta_max) {
      beta_min <- -1; beta_max <- 1
    }

    p <- plotParams()

    res <- as.numeric(input$resolution)
    beta_vals <- seq(beta_min, beta_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = p$mu_gamma, var = p$sigma2_gamma),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    post_vals <- gamma_beta_post(p$gamma_val, beta_vals, data, prior, log = FALSE, verbose = FALSE)
    df <- data.frame(beta = beta_vals, posterior = as.numeric(post_vals))

    # find beta modes
    beta_modes <- find_all_roots(beta_marg_post_drv,
                                 lower = beta_min, upper = beta_max,
                                 gamma = p$gamma_val, data = data, prior = prior,
                                 n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                                 eps_small = 1e-8)
    if (length(beta_modes) > 1) beta_modes <- beta_modes[-2]

    ggplot(df, aes(x = beta, y = posterior)) +
      geom_line(color = "darkorange", linewidth = 1.2) +
      geom_vline(xintercept = p$mu_beta, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = beta_modes, linetype = "dashed", color = "black") +
      labs(x = expression(beta), y = paste("Marginal Posterior at γ =", p$gamma_val)) +
      theme_minimal(base_size = 14)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
