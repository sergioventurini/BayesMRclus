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
