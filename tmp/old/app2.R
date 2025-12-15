# Functions
logpost_beta <- function(beta, gamma, prior, data, log = TRUE) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]
  
  tau2_j <- sigma2_Y + tau2
  
  loglik_Gammahat <- numeric(length(beta))
  for (i in 1:length(beta)) {
    h2_j <- (beta[i]^2)*psi2 + tau2_j
    loglik_Gammahat[i] <- -0.5*sum(log(h2_j) + (Gammahat_j - beta[i]*gamma)^2/h2_j)
    # loglik_Gammahat[i] <- sum(dnorm(Gammahat_j, mean = beta[i]*gamma, sd = sqrt(h2_j), log = TRUE))
  }
  logprior_beta <- -0.5*(beta - mu_beta)^2/sigma2_beta
  # logprior_beta <- dnorm(beta, mean = mu_beta, sd = sqrt(sigma2_beta), log = TRUE)

  # data.frame(loglik_Gammahat = loglik_Gammahat, logprior_beta = logprior_beta,
  #            lopost_beta = loglik_Gammahat + logprior_beta)
  res <- loglik_Gammahat + logprior_beta

  if (!log)
    res <- exp(res)

  res
}

logpost_beta_util <- function(beta, gamma, prior, data) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]
  
  tau2_j <- sigma2_Y + tau2
  
  loglik_Gammahat <- numeric(length(beta))
  for (i in 1:length(beta)) {
    h2_j <- (beta[i]^2)*psi2 + tau2_j
    loglik_Gammahat[i] <- -0.5*sum(log(h2_j) + (Gammahat_j - beta[i]*gamma)^2/h2_j)
    # loglik_Gammahat[i] <- sum(dnorm(Gammahat_j, mean = beta[i]*gamma, sd = sqrt(h2_j), log = TRUE))
  }
  logprior_beta <- -0.5*(beta - mu_beta)^2/sigma2_beta
  # logprior_beta <- dnorm(beta, mean = mu_beta, sd = sqrt(sigma2_beta), log = TRUE)

  data.frame(loglik_Gammahat = loglik_Gammahat, logprior_beta = logprior_beta)
}

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
      h2_j <- beta[b]^2*psi2 + tau2_j
      p_Gj <- log(h2_j) + (Gammahat_j - beta[b]*gamma[g])^2/h2_j
      
      res[g, b] <- -0.5*(sum(p_gj) + sum(p_Gj) + p_gamma + p_beta)
    }
  }
  res <- res[, , drop = TRUE]

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

    res[b] <- beta[b]^3*psi2^2*sum(1/h_j_beta^2) +
              beta[b]^2*(gamma*psi2*sum(Gammahat_j/h_j_beta^2)) +
              beta[b]*(gamma^2*sum(tau2_j/h_j_beta^2) - psi2*sum(Gammahat_j^2/h_j_beta^2) +
                       psi2*sum(tau2_j/h_j_beta^2) + 1/sigma2_beta) -
              (gamma*sum(Gammahat_j*tau2_j/h_j_beta^2) + mu_beta/sigma2_beta)
  }

  res
}

gamma_marg_post_mode <- function(beta, data, prior) {
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
  modes <- B_beta/A_beta

  res <- list(modes = modes,
              dens = gamma_beta_post(modes, beta, data, prior, log = FALSE, verbose = FALSE))

  res
}

find_all_roots <- function(f, ..., lower, upper, n = 1000,
                           tol_x = 1e-10, tol_f = 1e-12,
                           eps_small = 1e-8,
                           expand = FALSE,
                           step_out = 10,
                           max_steps = 100) {
  
  find_roots_in_interval <- function(a, b) {
    xs <- seq(a, b, length.out = n)
    vals <- sapply(xs, f, ...)
    
    brackets <- list()
    near_zero_intervals <- list()
    
    for (i in seq_len(n - 1)) {
      if (vals[i] == 0) {
        brackets[[length(brackets) + 1]] <- c(xs[i], xs[i])
      } else if (vals[i] * vals[i + 1] < 0) {
        brackets[[length(brackets) + 1]] <- c(xs[i], xs[i + 1])
      } else if (min(abs(vals[i]), abs(vals[i + 1])) < eps_small) {
        near_zero_intervals[[length(near_zero_intervals) + 1]] <- c(xs[i], xs[i + 1])
      }
    }
    
    refine_root <- function(a, b) {
      if (a == b) return(a)
      tryCatch({
        uniroot(f, ..., lower = a, upper = b, tol = tol_x)$root
      }, error = function(e) NA)
    }
    
    roots <- numeric(0)
    for (br in brackets) roots <- c(roots, refine_root(br[1], br[2]))
    for (br in near_zero_intervals) {
      mid <- mean(br)
      if (abs(f(mid, ...)) < tol_f * 100) {
        roots <- c(roots, refine_root(br[1], br[2]))
      }
    }
    
    roots <- roots[!is.na(roots)]
    sort(unique(roots))
  }
  
  # Initial search
  all_roots <- find_roots_in_interval(lower, upper)
  
  if (expand) {
    current_lower <- lower
    current_upper <- upper
    no_new_lower <- FALSE
    no_new_upper <- FALSE
    steps <- 0
    
    while ((!no_new_lower || !no_new_upper) && steps < max_steps) {
      steps <- steps + 1
      
      # Expand left
      if (!no_new_lower) {
        new_lower <- current_lower - step_out
        roots_left <- find_roots_in_interval(new_lower, current_lower)
        roots_left <- roots_left[roots_left < current_lower]
        if (length(roots_left) == 0) {
          no_new_lower <- TRUE
        } else {
          all_roots <- c(all_roots, roots_left)
          current_lower <- new_lower
        }
      }
      
      # Expand right
      if (!no_new_upper) {
        new_upper <- current_upper + step_out
        roots_right <- find_roots_in_interval(current_upper, new_upper)
        roots_right <- roots_right[roots_right > current_upper]
        if (length(roots_right) == 0) {
          no_new_upper <- TRUE
        } else {
          all_roots <- c(all_roots, roots_right)
          current_upper <- new_upper
        }
      }
    }
  }
  
  sort(unique(all_roots))
}

beta_marg_post_mode <- function(gamma, data, prior, beta_min = -20, beta_max = 20,
                                beta_step = 0.001, n = 1000, tol_x = 1e-10,
                                tol_f = 1e-12, eps_small = 1e-8, log = FALSE) {
  modes <- find_all_roots(beta_marg_post_drv,
                          lower = beta_min, upper = beta_max,
                          gamma = gamma, data = data, prior = prior,
                          n = n, tol_x = tol_x, tol_f = tol_f,
                          eps_small = eps_small, expand = TRUE,  # let the search automatically expand outside limits
                          step_out = 10, max_steps = 100)
  if (length(modes) > 1) modes <- modes[-2]

  res <- list(modes = modes,
              dens = gamma_beta_post(gamma, modes, data, prior, log = log, verbose = FALSE))

  res
}

bayesmr_noclus_optim <- function(data, prior, start = rep(0, 2), maxiter = 1000, tol = 1e-10,
                                 beta_min = -20, beta_max = 20, beta_step = 0.001,
                                 n = 1000, tol_x = 1e-10, tol_f = 1e-12, eps_small = 1e-8,
                                 verbose = TRUE) {
  gamma_chain <- beta_chain <- logpost <- NA

  # start iterations
  if (verbose) message("Running the optimization procedure...")
  
  gamma_chain <- start[1]  # gamma starting value
  beta_chain <- start[2]   # beta starting value
  logpost <- gamma_beta_post(gamma_chain[1], beta_chain[1], data, prior, log = TRUE, verbose = FALSE)
  if (verbose) message("  - iteration ", 1, " - gamma value: ", round(gamma_chain[1], digits = 6),
                       " - beta value: ", round(beta_chain[1], digits = 6),
                       " - log posterior value: ", format(logpost, nsamll = 8, scientific = FALSE))
  for (i in 2:maxiter) {
    # find gamma mode conditionally on beta
    gamma_optim <- gamma_marg_post_mode(beta_chain[i - 1], data, prior)
    gamma_chain <- c(gamma_chain, gamma_optim$modes)

    # find beta mode(s) conditionally on gamma
    beta_optim <- beta_marg_post_mode(gamma_chain[i], data, prior,
                                      beta_min = beta_min, beta_max = beta_max,
                                      beta_step = beta_step, n = n, tol_x = tol_x,
                                      tol_f = tol_f, eps_small = eps_small, log = TRUE)
    if (length(beta_optim$modes) == 1) {
      beta_chain <- c(beta_chain, beta_optim$modes)
    }
    else {
      if (identical(beta_optim$dens[1], beta_optim$dens[2]))
        stop("the joint posterior has two modes.")
      beta_chain <- c(beta_chain, beta_optim$modes[which.max(beta_optim$dens)])
    }

    logdens <- gamma_beta_post(gamma_chain[i], beta_chain[i], data, prior, log = TRUE, verbose = FALSE)
    logpost <- c(logpost, logdens)

    if (verbose) message("  - iteration ", i, " - gamma value: ", round(gamma_chain[i], digits = 6),
                         " - beta value: ", round(beta_chain[i], digits = 6),
                         " - log posterior value: ", format(logdens, nsamll = 8, scientific = FALSE))

    # check convergence
    if (max(abs(diff(gamma_chain[(i - 1):i])), abs(diff(beta_chain[(i - 1):i]))) < tol) {
      break
    }
  }

  # return results
  out <- list(gamma_best = gamma_chain[i], beta_best = beta_chain[i],
              gamma_chain = gamma_chain, beta_chain = beta_chain,
              logpost = logpost, iter = i)

  out
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
  titlePanel("Interactive (Unnormalized) Beta Marginal Posterior"),

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
    column(3, tags$div("$$\\beta \\text{ min}$$"), textInput("beta_min", NULL, value = "-1")),
    column(3, tags$div("$$\\beta \\text{ max}$$"), textInput("beta_max", NULL, value = "1"))
  ),

  fluidRow(column(12, tags$h4("Prior Parameters", style = "margin-top: 30px;"))),

  fluidRow(
    slider_cols(3, "mu_beta", "$$\\mu_\\beta$$", -5, 5, 0, step = 0.1),
    slider_cols(3, "sigma2_beta", "$$\\sigma^2_\\beta$$", 0, 10000, 0.1, step = 0.01),
    slider_cols(3, "psi2", "$$\\psi^2$$", 0, 10, 0.1, step = 0.01),
    slider_cols(3, "tau2", "$$\\tau^2$$", 0, 10, 0.1, step = 0.01)
  ),

  fluidRow(
    slider_cols(3, "gamma_for_beta", "$$\\gamma \\text{ value}$$", -20, 20, 0, step = 0.01)
  ),

  fluidRow(
    column(6, plotOutput("beta_plot", height = "300px")),
    column(6, plotOutput("beta_comp_plot", height = "300px"))
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
      mu_beta      = input$mu_beta,
      sigma2_beta  = input$sigma2_beta,
      psi2         = input$psi2,
      tau2         = input$tau2,
      gamma_val    = input$gamma_for_beta
    )
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

    res <- 300
    beta_vals <- seq(beta_min, beta_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = 0, var = 0.1),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    # post_vals <- gamma_beta_post(p$gamma_val, beta_vals, data, prior, log = FALSE, verbose = FALSE)
    post_vals <- logpost_beta(beta_vals, p$gamma_val, prior, data, log = FALSE)
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
      # geom_vline(xintercept = p$mu_beta, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = beta_modes, linetype = "dashed", color = "black") +
      labs(x = expression(beta), y = paste("Marginal Posterior at γ =", p$gamma_val)) +
      theme_minimal(base_size = 14)
  })

  output$beta_comp_plot <- renderPlot({
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

    res <- 300
    beta_vals <- seq(beta_min, beta_max, length.out = res)

    prior <- bayesmr_prior(
      gammaj = list(psi2 = p$psi2),
      Gammaj = list(tau2 = p$tau2),
      gamma  = list(mean = 0, var = 0.1),
      beta   = list(mean = p$mu_beta, var = p$sigma2_beta)
    )

    out <- logpost_beta_util(beta_vals, p$gamma_val, prior, data)
    df <- data.frame(beta = beta_vals, loglik = as.numeric(out[, 1]), logprior = as.numeric(out[, 2]))
    df$logpost <- df$loglik + df$logprior
    df_long <- df |>
      tidyr::pivot_longer(cols = c("loglik", "logprior", "logpost"),
                   names_to = "component",
                   values_to = "value")

    # find beta modes
    beta_modes <- find_all_roots(beta_marg_post_drv,
                                 lower = beta_min, upper = beta_max,
                                 gamma = p$gamma_val, data = data, prior = prior,
                                 n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                                 eps_small = 1e-8)
    if (length(beta_modes) > 1) beta_modes <- beta_modes[-2]

    ggplot(df_long, aes(x = beta, y = value, color = component)) +
      geom_line(linewidth = 1.2) +
      # geom_vline(xintercept = p$mu_beta, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = beta_modes, linetype = "dashed", color = "black") +
      labs(x = expression(beta), 
           y = paste("Marginal Log Curves at γ =", p$gamma_val),
           color = "Curve") + 
      theme_minimal(base_size = 14)
  })
}

# Run the app
shinyApp(ui = ui, server = server)
