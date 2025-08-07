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

# UI
ui <- fluidPage(
  titlePanel("Interactive (Unnormalized) Joint Posterior"),

  withMathJax(),

  fluidRow(column(12, tags$h4("Axis Ranges", style = "margin-top: 20px;"))),
  fluidRow(
    column(3, tags$div("$$\\gamma \\text{ min}$$"), textInput("gamma_min", NULL, value = "-0.05")),
    column(3, tags$div("$$\\gamma \\text{ max}$$"), textInput("gamma_max", NULL, value = "0.05")),
    column(3, tags$div("$$\\beta \\text{ min}$$"), textInput("beta_min", NULL, value = "-1")),
    column(3, tags$div("$$\\beta \\text{ max}$$"), textInput("beta_max", NULL, value = "1"))
  ),

  fluidRow(column(12, tags$h4("Prior Parameters", style = "margin-top: 30px;"))),
  fluidRow(
    column(4, tags$div("$$\\mu_\\gamma$$"), sliderInput("mu_gamma", NULL, -5, 5, 0, step = 1e-2)),
    column(4, tags$div("$$\\mu_\\beta$$"), sliderInput("mu_beta", NULL, -5, 5, 0, step = 1e-2)),
    column(4, tags$div("$$\\psi^2$$"), sliderInput("psi2", NULL, 0.0001, 10, 0.1, step = 1e-3))
  ),
  fluidRow(
    column(4, tags$div("$$\\sigma^2_\\gamma$$"), sliderInput("sigma2_gamma", NULL, 0.0001, 10, 0.1, step = 1e-3)),
    column(4, tags$div("$$\\sigma^2_\\beta$$"), sliderInput("sigma2_beta", NULL, 0.0001, 10, 0.1, step = 1e-3)),
    column(4, tags$div("$$\\tau^2$$"), sliderInput("tau2", NULL, 0.0001, 10, 0.1, step = 1e-3))
  ),

  fluidRow(column(12, plotOutput("plot", height = "700px")))
)

# Server
server <- function(input, output, session) {
  output$plot <- renderPlot({
    # Validate numeric ranges
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

    # Parameters
    mu_gamma     <- input$mu_gamma
    sigma2_gamma <- input$sigma2_gamma
    mu_beta      <- input$mu_beta
    sigma2_beta  <- input$sigma2_beta
    psi2         <- input$psi2
    tau2         <- input$tau2

    # Grid
    gamma_vals <- seq(gamma_min, gamma_max, length.out = 150)
    beta_vals  <- seq(beta_min, beta_max, length.out = 150)

    # Posterior
    prior <- bayesmr_prior(
      gammaj = list(psi2 = psi2),
      Gammaj = list(tau2 = tau2),
      gamma  = list(mean = mu_gamma, var = sigma2_gamma),
      beta   = list(mean = mu_beta, var = sigma2_beta)
    )
    post_vals <- gamma_beta_post(gamma_vals, beta_vals, data, prior, log = FALSE, verbose = FALSE)

    # Convert to data frame for ggplot
    df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
    df_plot$posterior <- as.vector(post_vals)

    # Plot
    ggplot(df_plot, aes(x = gamma, y = beta, z = posterior)) +
      geom_contour_filled(bins = 20) +
      scale_fill_viridis_d(option = "C") +
      geom_vline(xintercept = mu_gamma, linetype = "dashed", color = "gray") +
      geom_hline(yintercept = mu_beta, linetype = "dashed", color = "gray") +
      labs(
        x = expression(gamma),
        y = expression(beta),
        fill = "Posterior"
      ) +
      coord_cartesian(xlim = c(gamma_min, gamma_max),
                      ylim = c(beta_min, beta_max)) +
      theme_minimal(base_size = 14)
  })
}

# Run app
shinyApp(ui = ui, server = server)
