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
profile_loglik <- function(beta, psi2, tau2, data) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  # stacked vector
  v <- c(gammahat_j, Gammahat_j)

  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  a <- 1/psi2_j
  b <- 1/(beta^2*psi2 + tau2_j)

  # diagonal part
  D <- c(a, b)

  # rank-one part
  u <- c(a, beta*b)
  s <- sum(a) + beta^2*sum(b)

  # quadratic form: v^T (D - uu^T/s) v
  quad <- sum(D*v^2) - (sum(u*v)^2)/s

  res <- -0.5*(quad - sum(log(b)))

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
profile_loglik_all <- function(beta, psi2, tau2, data) {
  beta_len <- length(beta)
  psi2_len <- length(psi2)
  tau2_len <- length(tau2)

  res <- array(NA, dim = c(beta_len, psi2_len, tau2_len))
  for (b in 1:beta_len) {
    for (p in 1:psi2_len) {
      for (t in 1:tau2_len) {
        res[b, p, t] <- profile_loglik(beta[b], psi2[p], tau2[t], data)
      }
    }
  }

  drop(res)
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
fit_profile_loglik <- function(data, init = c(beta = 0, logpsi2 = -2, logtau2 = -2)) {
  nll <- function(par) {
    beta <- par[1]
    psi2 <- exp(par[2])  # ensure positivity
    tau2 <- exp(par[3])  # ensure positivity
    -profile_loglik(beta, psi2, tau2, data)
  }
  
  opt <- optim(
    par = init,
    fn = nll,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 1000)
  )
  
  # Extract estimates
  beta_hat <- opt$par[1]
  psi2_hat <- exp(opt$par[2])
  tau2_hat <- exp(opt$par[3])
  loglik   <- -opt$value
  
  # Compute covariance matrix (inverse Hessian)
  cov_par <- tryCatch(
    solve(opt$hessian),
    error = function(e) matrix(NA, 3, 3)
  )
  
  # Standard errors
  se_beta <- sqrt(cov_par[1, 1])
  
  # Delta method for psi2: Var(psi2) ≈ (d psi2 / d logpsi2)^2 Var(logpsi2)
  se_psi2 <- psi2_hat * sqrt(cov_par[2, 2])

  # Delta method for tau2: Var(tau2) ≈ (d tau2 / d logtau2)^2 Var(logtau2)
  se_tau2 <- tau2_hat * sqrt(cov_par[3, 3])

  list(beta = as.numeric(beta_hat),
       psi2 = as.numeric(psi2_hat),
       tau2 = as.numeric(tau2_hat),
       se_beta = as.numeric(se_beta),
       se_psi2 = as.numeric(se_psi2),
       se_tau2 = as.numeric(se_tau2),
       cov_par = cov_par,
       loglik = loglik,
       opt = opt)
}
