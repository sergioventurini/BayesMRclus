#' Function to implement the Metropolis algorithm for an arbitrary posterior probability distribution
#'
#' \code{metropolis()} implements a general Metropolis MCMC algorithm that can be applied to any
#' posterior probability provided by the user.
#'
#' @param logpost A function defining the logarithm of the posterior density.
#' @param current A numeric value representing the starting value.
#' @param proposal A length-one character vector with the name of the proposal
#'   distribution; currently, accepted values are "unif" or "norm".
#' @param C A numeric value corresponding to the neighborhood where one looks
#'   for a proposal value; for a uniform proposal is the distance from the
#'   current value; for a normal proposal it is the standard deviation.
#' @param iter Integer value specifying the total the number of iterations
#'   of the algorithm.
#' @param ... Additional arguments to be passed to the \code{logpost}() function.
#' @return A list with two components, \code{S} is a vector of the simulated draws,
#' and \code{accept_rate}, which gives the acceptance rate of the algorithm.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{TODO}()} for computing ...
#' @examples
#' lpost <- function(theta, s) {
#'   dcauchy(theta, s$loc, s$scale, log = TRUE) +
#'   dnorm(s$ybar, theta, s$se, log = TRUE)
#' }
#' 
#' data <- read.table("./data/buffalo_jan.txt", header = TRUE)
#' s <- list(loc = 10, scale = 2,
#'           ybar = mean(data$JAN),
#'           se = sd(data$JAN) / sqrt(nrow(data)))
#'
#' set.seed(101)
#' out <- metropolis(lpost, 5, 20, 10000, s)
#'
#' @export
metropolis <- function(logpost, current, proposal = "unif", C, iter, ...) {
  S <- rep(0, iter)
  n_accept <- 0
  for (j in 1:iter) {
    # uniform proposal distribution centred around the current value
    if (proposal == "unif") {
      candidate <- runif(1, min = current - C, max = current + C)
    }
    else if (proposal == "norm") {
      candidate <- rnorm(1, mean = current, sd = C)
    }
    else {
      stop("the specified proposal distribution is not available")
    }
    prob <- exp(logpost(candidate, ...) - logpost(current, ...))
    accept <- ifelse(runif(1) < prob, 1, 0)
    current <- ifelse(accept == 1, candidate, current)
    S[j] <- current
  }
  n_accept <- sum(accept)

  res <- list(S = S, accept_rate = n_accept / iter)

  return(res)
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
logpost_gamma <- function(g, par, hpar, data) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  beta <- par[["beta"]]

  mu_gamma <- hpar[["mu_gamma"]]
  sigma2_gamma <- hpar[["sigma2_gamma"]]
  psi2 <- hpar[["psi2"]]
  tau2 <- hpar[["tau2"]]

  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2
  h2_j <- (beta^2)*psi2 + tau2_j

  A_beta <- sum(1/psi2_j) + (beta^2)*sum(1/h2_j) + 1/sigma2_gamma
  B_beta <- sum(gammahat_j/psi2_j) + beta*sum(Gammahat_j/h2_j) + mu_gamma/sigma2_gamma

  res <- dnorm(g, B_beta/A_beta, sqrt(1/A_beta), log = TRUE)

  return(res)
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
#' hpar <- list(mu_beta = 0,
#'              sigma2_beta = 0.1,
#'              psi2 = 0.3,
#'              tau2 = 0.4)
#'
#' b <- seq(-0.8, 0.8, length.out = 100)
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#' par <- list(gamma = 0.5)
#' res <- logpost_beta(b, par, hpar, data)
#' summary(res)
#'
#' @export
logpost_beta <- function(b, par, hpar, data) {
  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  gamma <- par[["gamma"]]

  mu_beta <- hpar[["mu_beta"]]
  sigma2_beta <- hpar[["sigma2_beta"]]
  psi2 <- hpar[["psi2"]]
  tau2 <- hpar[["tau2"]]
  
  tau2_j <- sigma2_Y + tau2
  
  loglik_Gammahat <- numeric(length(b))
  for (i in 1:length(b)) {
    h2_j <- (b[i]^2)*psi2 + tau2_j
    loglik_Gammahat[i] <- -0.5*(sum(log(h2_j) + (Gammahat_j - b[i]*gamma)^2/h2_j))
  }
  logprior_beta <- -0.5*(b - mu_beta)^2/sigma2_beta
  res <- loglik_Gammahat + logprior_beta

  return(res)
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
#' hpar <- list(mu_beta = 0,
#'              sigma2_beta = 0.1,
#'              psi2 = 0.3,
#'              tau2 = 0.4)
#'
#' b <- seq(-0.8, 0.8, length.out = 100)
#' data <- bmi_sbp[, c("beta.exposure",
#'                     "beta.outcome",
#'                     "se.exposure",
#'                     "se.outcome")]
#' par <- list(gamma = 0.5)
#' res <- logpost_beta(b, par, hpar, data)
#' summary(res)
#'
#' @export
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
  
  loglik_Gammahat <- numeric(length(b))
  for (i in 1:length(beta)) {
    h2_j <- (beta[i]^2)*psi2 + tau2_j
    loglik_Gammahat[i] <- -0.5*(sum(log(h2_j)) + sum((Gammahat_j - beta[i]*gamma)^2/h2_j))
  }
  logprior_beta <- -0.5*(beta - mu_beta)^2/sigma2_beta

  data.frame(loglik_Gammahat = loglik_Gammahat, logprior_beta = logprior_beta)
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
mcmc_bayesmr <- function(data, hpar, iter, start, tune, proposal = "norm") {
  if (is.null(proposal)) {
    proposal <- "unif"
  }

  gammahat_j <- data[, 1]  # SNP-Exposure effect
  Gammahat_j <- data[, 2]  # SNP-Outcome effect
  sigma2_X <- data[, 3]^2  # SNP-Exposure effect variance
  sigma2_Y <- data[, 4]^2  # SNP-Outcome effect variance

  mu_gamma <- hpar[["mu_gamma"]]
  sigma2_gamma <- hpar[["sigma2_gamma"]]
  mu_beta <- hpar[["mu_beta"]]
  sigma2_beta <- hpar[["sigma2_beta"]]
  psi2 <- hpar[["psi2"]]
  tau2 <- hpar[["tau2"]]

  psi2_j <- sigma2_X + psi2
  tau2_j <- sigma2_Y + tau2

  gamma_p <- start[1]
  beta_p <- start[2]

  draws <- data.frame(gamma = numeric(0), beta = numeric(0))
  C <- tune
  accept_rate <- 0
  for (m in 1:iter) {
    # gamma update using its full conditional
    h2_j <- (beta_p^2)*psi2 + tau2_j
    A_beta <- sum(1/psi2_j) + (beta_p^2)*sum(1/h2_j) + 1/sigma2_gamma
    B_beta <- sum(gammahat_j/psi2_j) + beta_p*sum(Gammahat_j/h2_j) + mu_gamma/sigma2_gamma
    gamma_p <- rnorm(1, B_beta/A_beta, sqrt(1/A_beta))
    draws[m, 1] <- gamma_p

    # beta update using a M-H step
    par <- list(gamma = gamma_p)
    hpar <- list(mu_beta = mu_beta, sigma2_beta = sigma2_beta, psi2 = psi2, tau2 = tau2)
    mh <- metropolis(logpost = logpost_beta, current = beta_p, proposal = proposal,
                     C = C, iter = 1, par, hpar, data)
    beta_p <- mh$S
    draws[m, 2] <- beta_p
    accept_rate <- accept_rate + mh$accept_rate

    # if ((m %% 500) == 0) print(paste0("Simulation ", m, " of ", format(iter)))
  }

  res <- list(draws = draws, accept_rate = accept_rate/iter)

  return(res)
}
