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
mcmc_bayesmr_profile <- function(data, prior, iter, start, tune, verbose = TRUE) {
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

  current_beta <- start[1]
  current_prec_psi2 <- 1/start[2]
  current_prec_tau2 <- 1/start[3]

  draws <- data.frame(beta = numeric(0), psi2 = numeric(0), tau2 = numeric(0))
  accept_beta <- accept_psi2 <- accept_tau2 <- 0
  accept_rate <- numeric(3)
  for (m in 1:iter) {
    # beta update using a M-H step
    beta_p <- rnorm(1, mean = current_beta, sd = tune[1])
    prob_beta <- exp(logpost_profile(beta_p, current_prec_psi2, current_prec_tau2, prior, data) -
                     logpost_profile(current_beta, current_prec_psi2, current_prec_tau2, prior, data))
    accept_beta <- ifelse(runif(1) < prob_beta, 1, 0)
    current_beta <- ifelse(accept_beta, beta_p, current_beta)
    draws[m, 1] <- current_beta
    accept_rate[1] <- accept_rate[1] + accept_beta

    # prec_psi2 update using a M-H step
    prec_psi2_p <- runif(1, max(c(1e-1, current_prec_psi2 - tune[2])), min(c(1e3, current_prec_psi2 + tune[2])))
    prob_prec_psi2 <- exp(logpost_profile(current_beta, prec_psi2_p, current_prec_tau2, prior, data) -
                          logpost_profile(current_beta, current_prec_psi2, current_prec_tau2, prior, data))
    accept_prec_psi2 <- ifelse(runif(1) < prob_prec_psi2, 1, 0)
    current_prec_psi2 <- ifelse(accept_prec_psi2, prec_psi2_p, current_prec_psi2)
    draws[m, 2] <- 1/current_prec_psi2
    accept_rate[2] <- accept_rate[2] + accept_prec_psi2

    # prec_tau2 update using a M-H step
    prec_tau2_p <- runif(1, max(c(1e-1, current_prec_tau2 - tune[3])), min(c(1e8, current_prec_tau2 + tune[3])))
    prob_prec_tau2 <- exp(logpost_profile(current_beta, current_prec_psi2, prec_tau2_p, prior, data) -
                          logpost_profile(current_beta, current_prec_psi2, current_prec_tau2, prior, data))
    accept_prec_tau2 <- ifelse(runif(1) < prob_prec_tau2, 1, 0)
    current_prec_tau2 <- ifelse(accept_prec_tau2, prec_tau2_p, current_prec_tau2)
    draws[m, 3] <- 1/current_prec_tau2
    accept_rate[3] <- accept_rate[3] + accept_prec_tau2

    if (verbose)
      if ((m %% 500) == 0) print(paste0("Simulation ", m, " of ", format(iter)))
  }

  res <- list(draws = draws, accept_rate = accept_rate/iter)

  return(res)
}
