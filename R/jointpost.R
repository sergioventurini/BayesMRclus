#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param data Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param prior A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value.
#' @param iter Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param start A numeric vector providing the parameter starting values.
#' @param tune A numeric vector providing the Metropolis-Hastings tuning parameters.
#' @param proposal A length-one character vector with the name of the proposal
#'   distribution; currently, accepted values are "unif" or "norm".
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' data(bmi_sbp)
#' 
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#'
#' beta_z <- 0.3
#' beta_se_z <- 0.2
#' tau2_z <- 1e-3
#' hpar <- list(mu_gamma = mean(data[, 1]),
#'              sigma2_gamma = var(data[, 1])/nrow(data),
#'              mu_beta = beta_z,
#'              sigma2_beta = beta_se_z^2,
#'              psi2 = mean(data[, 3]^2),
#'              tau2 = tau2_z)
#'
#' iter <- 1e4
#' start <- rnorm(2)
#' tune <- 1.5
#' res <- mcmc_bayesmr(data, hpar, iter, start, tune)
#' summary(res$draws)
#' res$accept_rate
#'
#' @export
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
  res <- res[, , drop = TRUE]

  if (!log)
    res <- exp(res)

  res
}

#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param data Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param prior A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value.
#' @param iter Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param start A numeric vector providing the parameter starting values.
#' @param tune A numeric vector providing the Metropolis-Hastings tuning parameters.
#' @param proposal A length-one character vector with the name of the proposal
#'   distribution; currently, accepted values are "unif" or "norm".
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' data(bmi_sbp)
#' 
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#'
#' beta_z <- 0.3
#' beta_se_z <- 0.2
#' tau2_z <- 1e-3
#' hpar <- list(mu_gamma = mean(data[, 1]),
#'              sigma2_gamma = var(data[, 1])/nrow(data),
#'              mu_beta = beta_z,
#'              sigma2_beta = beta_se_z^2,
#'              psi2 = mean(data[, 3]^2),
#'              tau2 = tau2_z)
#'
#' iter <- 1e4
#' start <- rnorm(2)
#' tune <- 1.5
#' res <- mcmc_bayesmr(data, hpar, iter, start, tune)
#' summary(res$draws)
#' res$accept_rate
#'
#' @export
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

#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param data Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param prior A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value.
#' @param iter Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param start A numeric vector providing the parameter starting values.
#' @param tune A numeric vector providing the Metropolis-Hastings tuning parameters.
#' @param proposal A length-one character vector with the name of the proposal
#'   distribution; currently, accepted values are "unif" or "norm".
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' data(bmi_sbp)
#' 
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#'
#' beta_z <- 0.3
#' beta_se_z <- 0.2
#' tau2_z <- 1e-3
#' hpar <- list(mu_gamma = mean(data[, 1]),
#'              sigma2_gamma = var(data[, 1])/nrow(data),
#'              mu_beta = beta_z,
#'              sigma2_beta = beta_se_z^2,
#'              psi2 = mean(data[, 3]^2),
#'              tau2 = tau2_z)
#'
#' iter <- 1e4
#' start <- rnorm(2)
#' tune <- 1.5
#' res <- mcmc_bayesmr(data, hpar, iter, start, tune)
#' summary(res$draws)
#' res$accept_rate
#'
#' @export
beta_marg_post_mode <- function(gamma, data, prior, beta_min = -20, beta_max = 20,
                                beta_step = 0.001, n = 1000, tol_x = 1e-10,
                                tol_f = 1e-12, eps_small = 1e-8, log = FALSE) {
  modes <- find_all_roots(beta_marg_post_drv,
                          lower = beta_min, upper = beta_max,
                          gamma = gamma, data = data, prior = prior,
                          n = n, tol_x = tol_x, tol_f = tol_f,
                          eps_small = eps_small, expand = TRUE,  # let the search automatically expand outside limits
                          step_out = 10, max_steps = 100)
  if (length(modes) > 1) {
    modes <- modes[-2]
  }

  res <- list(modes = modes,
              dens = gamma_beta_post(gamma, modes, data, prior, log = log, verbose = FALSE))

  res
}

#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param data Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param hpar A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value.
#' @param iter Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param start A numeric vector providing the parameter starting values.
#' @param tune A numeric vector providing the Metropolis-Hastings tuning parameters.
#' @param proposal A length-one character vector with the name of the proposal
#'   distribution; currently, accepted values are "unif" or "norm".
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' data(bmi_sbp)
#' 
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#'
#' beta_z <- 0.3
#' beta_se_z <- 0.2
#' tau2_z <- 1e-3
#' hpar <- list(mu_gamma = mean(data[, 1]),
#'              sigma2_gamma = var(data[, 1])/nrow(data),
#'              mu_beta = beta_z,
#'              sigma2_beta = beta_se_z^2,
#'              psi2 = mean(data[, 3]^2),
#'              tau2 = tau2_z)
#'
#' iter <- 1e4
#' start <- rnorm(2)
#' tune <- 1.5
#' res <- mcmc_bayesmr(data, hpar, iter, start, tune)
#' summary(res$draws)
#' res$accept_rate
#'
#' @export
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

#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param data Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param prior A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value.
#' @param iter Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param start A numeric vector providing the parameter starting values.
#' @param tune A numeric vector providing the Metropolis-Hastings tuning parameters.
#' @param proposal A length-one character vector with the name of the proposal
#'   distribution; currently, accepted values are "unif" or "norm".
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' data(bmi_sbp)
#' 
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#'
#' beta_z <- 0.3
#' beta_se_z <- 0.2
#' tau2_z <- 1e-3
#' hpar <- list(mu_gamma = mean(data[, 1]),
#'              sigma2_gamma = var(data[, 1])/nrow(data),
#'              mu_beta = beta_z,
#'              sigma2_beta = beta_se_z^2,
#'              psi2 = mean(data[, 3]^2),
#'              tau2 = tau2_z)
#'
#' iter <- 1e4
#' start <- rnorm(2)
#' tune <- 1.5
#' res <- mcmc_bayesmr(data, hpar, iter, start, tune)
#' summary(res$draws)
#' res$accept_rate
#'
#' @export
beta_mode_dist <- function(gamma, data, prior,
                           tau2_min = 0, tau2_max = 5, tau2_step = 0.001,
                           beta_min = -20, beta_max = 20, beta_step = 0.001,
                           n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                           eps_small = 1e-8, verbose = TRUE) {
  tau2_seq <- seq(tau2_min, tau2_max, tau2_step)

  beta_modes_all <- list()
  beta_dist_all <- numeric(0)
  t <- 1
  beta_dist <- 1e10
  while (beta_dist > 0) {
    prior$Gammaj$tau2 <- tau2_seq[t]

    beta_modes <- beta_marg_post_mode(gamma, data, prior,
                                      beta_min = beta_min, beta_max = beta_max,
                                      beta_step = beta_step,
                                      n = n, tol_x = tol_x, tol_f = tol_f,
                                      eps_small = eps_small)
    
    roots <- beta_modes$modes
    beta_dist <- ifelse(length(roots) > 1, diff(roots), 0)
    beta_modes_all[[t]] <- beta_modes
    beta_dist_all <- c(beta_dist_all, beta_dist)
    
    if (verbose)
      message("tau2 value: ", format(tau2_seq[t], nsmall = -log10(tau2_step), scientific = FALSE),
              " - modes distance = ", round(beta_dist, digits = 6), " - iter = ", t)
    
    t <- t + 1
  }

  res <- list(tau2_star = tau2_seq[t],
              iter = t - 1,
              modes = beta_modes_all,
              dist = beta_dist_all)

  res
}
