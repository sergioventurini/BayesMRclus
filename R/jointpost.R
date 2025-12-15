#' Log-likelihood for BayesMR models.
#'
#' \code{bayesmr_logLik()} computes the log-likelihood value for a BayesMR model.
#'
#' @param gamma A length-one numeric vector containing the gamma parameter value
#' @param beta A length-one numeric vector containing the beta parameter value
#' @param data An object of class \code{\link{bayesmr_data}}).
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{bayesmr_prior}()} for more details.
#'
#' @return A numeric matrix of the log-likelihood values.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso \code{\link{bayesmr}()}.
#'
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#'
#' @export
bayesmr_logLik <- function(gamma, beta, data, prior, log = TRUE) {
  n <- nrow(data)

  gammaj_hat <- data[, 1]
  Gammaj_hat <- data[, 2]
  sigmaj_X <- data[, 3]
  sigmaj_Y <- data[, 4]

  # recover prior hyperparameters
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  psi2_j <- sigmaj_X^2 + psi2
  tau2_j <- sigmaj_Y^2 + tau2
  
  gamma_len <- length(gamma)
  beta_len <- length(beta)
  res <- matrix(NA, nrow = gamma_len, ncol = beta_len)
  ll <- numeric(n)

  for (g in 1:gamma_len) {
    for (b in 1:beta_len) {
      for (i in 1:n) {
        ll[i] <- dbivnorm(gammaj_hat[i], Gammaj_hat[i],
                          gamma[g], beta[b]*gamma[g],
                          psi2_j[i], beta[b]^2*psi2 + tau2_j[i], beta[b]*psi2,
                          log = TRUE)
      }
      res[g, b] <- sum(ll)
    }
  }
  res <- res[, , drop = TRUE]

  if (!log) res <- exp(res)

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
  # unpack data
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  # unpack prior
  mu_gamma <- prior[["gamma"]][["mean"]]
  sigma2_gamma <- prior[["gamma"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  # precompute per-SNP constants
  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  a_j <- (beta^2)*psi2 + tau2_j
  v_j <- a_j*psi2_j - beta*psi2^2

  A_beta <- sum(tau2_j/v_j + (beta^2)*sigma2_X/v_j) + 1/sigma2_gamma
  B_beta <- sum(tau2_j*gammahat_j/v_j + beta*sigma2_X*Gammahat_j/v_j) + mu_gamma/sigma2_gamma
  modes <- B_beta/A_beta

  res <- list(modes = modes,
              dens = gamma_beta_post(modes, beta, data, prior, log = FALSE))

  res
}

#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param g A function defining the logarithm of the posterior density.
#' @param par A numeric value representing the starting value.
#' @param hpar A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value.
#' @param data Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' data(bmi_sbp)
#' 
#' hpar <- list(mu_gamma = 0,
#'              sigma2_gamma = 1,
#'              psi2 = 3,
#'              tau2 = 4)
#'
#' g <- seq(-0.2, 0.2, length.out = 100)
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#' par <- list(beta = 0.3)
#' res <- logpost_gamma(g, par, hpar, data)
#' summary(res)
#'
#' @export
logpost_gamma <- function(gamma, beta, prior, data) {
  # unpack data
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  # unpack prior
  mu_gamma <- prior[["gamma"]][["mean"]]
  sigma2_gamma <- prior[["gamma"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  # precompute per-SNP constants
  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  a_j <- (beta^2)*psi2 + tau2_j
  v_j <- a_j*psi2_j - beta*psi2^2

  A_beta <- sum(tau2_j/v_j + (beta^2)*sigma2_X/v_j) + 1/sigma2_gamma
  B_beta <- sum(tau2_j*gammahat_j/v_j + beta*sigma2_X*Gammahat_j/v_j) + mu_gamma/sigma2_gamma

  res <- dnorm(gamma, B_beta/A_beta, sqrt(1/A_beta), log = TRUE)

  return(res)
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
