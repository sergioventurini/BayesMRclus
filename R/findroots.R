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
