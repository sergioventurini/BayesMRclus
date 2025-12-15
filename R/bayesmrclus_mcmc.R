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
metropolis <- function(logpost, current, proposal = c("unif", "norm"), C, iter, ...) {
  proposal <- match.arg(proposal)

  S <- numeric(iter)
  n_accept <- 0L

  # evaluate logpost(current) once per iteration (updated when accepted)
  logpost_current <- logpost(current, ...)
  for (j in seq_len(iter)) {
    # draw candidate
    candidate <- switch(proposal,
                        unif = runif(1, min = current - C, max = current + C),
                        norm = rnorm(1, mean = current, sd = C))

    # compute log ratio once
    logpost_candidate <- logpost(candidate, ...)
    log_ratio <- logpost_candidate - logpost_current

    # numerical stability & handle non-finite values:
    if (!is.finite(log_ratio)) {
      accept_prob <- 0
    } else {
      # acceptance probability = min(1, exp(log_ratio))
      # avoid calling exp on large positive log_ratio by using pmin
      accept_prob <- if (log_ratio <= 0) exp(log_ratio) else 1
    }

    # accept / reject
    if (runif(1) < accept_prob) {
      current <- candidate
      logpost_current <- logpost_candidate
      n_accept <- n_accept + 1L
    }

    S[j] <- current
  }

  list(S = S, accept_rate = n_accept / iter)
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
mcmc_bayesmr <- function(data, prior, iter, start, tune,
  proposal = c("unif", "norm"), verbose = TRUE) {
  proposal <- match.arg(proposal)

  # unpack data columns (assume they are in the order: gammahat, Gammahat, seX, seY)
  gammahat_j <- data[, 1]
  Gammahat_j <- data[, 2]
  seX_j <- data[, 3]
  seY_j <- data[, 4]

  sigma2_X <- seX_j^2
  sigma2_Y <- seY_j^2

  # unpack prior
  mu_gamma <- prior[["gamma"]][["mean"]]
  sigma2_gamma <- prior[["gamma"]][["var"]]
  mu_beta <- prior[["beta"]][["mean"]]
  sigma2_beta <- prior[["beta"]][["var"]]
  psi2 <- prior[["gammaj"]][["psi2"]]
  tau2 <- prior[["Gammaj"]][["tau2"]]

  # precompute per-SNP constants
  psi2_j <- sigma2_X + psi2    # vector
  tau2_j <- sigma2_Y + tau2    # vector

  # starting values
  gamma_p <- start[1]
  beta_p  <- start[2]

  # pre-allocate draws as matrix (faster than data.frame in loop)
  draws_mat <- matrix(NA_real_, nrow = iter, ncol = 2)
  colnames(draws_mat) <- c("gamma", "beta")

  C <- tune
  acc_sum <- 0      # will accumulate acceptance rates from metropolis (which returns proportion)

  for (m in seq_len(iter)) {
    ## --- gamma update: direct sample from conditional (as in original) ---
    # these use per-SNP vectors and current beta_p
    a_j <- (beta_p^2) * psi2 + tau2_j      # vector
    c_beta <- beta_p * psi2                # scalar
    v_j <- a_j * psi2_j - c_beta^2        # vector

    # avoid divide-by-zero / non-finite v_j (defensive)
    if (any(!is.finite(v_j)) || any(v_j <= 0)) {
      stop("non-finite or non-positive v_j encountered in gamma update.")
    }

    A_beta <- sum(tau2_j / v_j) + (beta_p^2) * sum(sigma2_X / v_j) + 1 / sigma2_gamma
    B_beta <- sum(tau2_j * gammahat_j / v_j) + beta_p * sum(sigma2_X * Gammahat_j / v_j) + mu_gamma / sigma2_gamma

    gamma_p <- rnorm(1, mean = B_beta / A_beta, sd = sqrt(1 / A_beta))
    draws_mat[m, 1] <- gamma_p

    # beta update using a M-H step
    ## --- beta update: single-step Metropolis (uses metropolis_opt with iter = 1) ---
    mh <- metropolis(logpost = logpost_beta,    # same external function name as original
                     current = beta_p,
                     proposal = proposal,
                     C = C,
                     iter = 1,
                     gamma = gamma_p,
                     prior = prior,
                     data = data)

    beta_p <- as.numeric(mh$S[1])
    draws_mat[m, 2] <- beta_p

    acc_sum <- acc_sum + mh$accept_rate

    if (verbose && (m %% 500L) == 0L) {
      message(sprintf("Simulation %d of %d", m, iter))
    }
  }

  list(draws = draws_mat, accept_rate = acc_sum / iter)
}
