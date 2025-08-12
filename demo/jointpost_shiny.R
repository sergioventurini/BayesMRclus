if (!requireNamespace("BayesMRclus", quietly = TRUE)) {
  remotes::install_github("sergioventurini/BayesMRclus")
}
library(BayesMRclus)

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
      choices = c("Low" = 50, "Medium" = 100, "High" = 150), selected = 50)),
    column(4,
        tags$div(
          style = "border: 1px solid #ddd; padding: 10px; background-color: #f9f9f9; border-radius: 5px;",
          "Note: the coloured dot and solid lines represent the joint posterior global maximum."))
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
    # beta_modes <- find_all_roots(beta_marg_post_drv,
    #                              lower = beta_min, upper = beta_max,
    #                              gamma = p$gamma_val, data = data, prior = prior,
    #                              n = 1000, tol_x = 1e-10, tol_f = 1e-12,
    #                              eps_small = 1e-8)
    # if (length(beta_modes) > 1) beta_modes <- beta_modes[-2]

    # find gamma mode
    # gamma_mode <- gamma_marg_post_mode(p$beta_val, data, prior)

    glob_max <- bayesmr_noclus_optim(data, prior, start = rep(0, 2), maxiter = 1000, tol = 1e-10,
                                     beta_min = beta_min, beta_max = beta_max, beta_step = 0.001,
                                     n = 1000, tol_x = 1e-10, tol_f = 1e-12, eps_small = 1e-8,
                                     verbose = FALSE)

    ggplot(df_plot, aes(x = gamma, y = beta, z = posterior)) +
      geom_contour_filled(bins = 20) +
      scale_fill_viridis_d(option = "C") +
      # geom_vline(xintercept = p$mu_gamma, linetype = "dashed", color = "gray") +
      # geom_vline(xintercept = gamma_mode$modes, linetype = "dashed", color = "black") +
      # geom_hline(yintercept = p$mu_beta, linetype = "dashed", color = "gray") +
      # geom_hline(yintercept = beta_modes, linetype = "dashed", color = "black") +
      geom_vline(xintercept = glob_max$gamma_best, linetype = "solid", color = "#21908C") +
      geom_hline(yintercept = glob_max$beta_best, linetype = "solid", color = "#21908C") +
      labs(x = expression(gamma), y = expression(beta), fill = "Posterior") +
      coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
      geom_point(aes(x = glob_max$gamma_best, y = glob_max$beta_best), color = "#21908C", size = 3) +
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
    gamma_mode <- gamma_marg_post_mode(p$beta_val, data, prior)

    ggplot(df, aes(x = gamma, y = posterior)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_vline(xintercept = p$mu_gamma, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = gamma_mode$modes, linetype = "dashed", color = "black") +
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
