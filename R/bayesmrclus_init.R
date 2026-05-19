#' Generate Starting Values for BayesMR Model Fitting
#'
#' \code{bayesmr_init()} draws or validates starting values for the MCMC
#' sampler before fitting a BayesMR model.
#'
#' @param data An object of class \code{\link{bayesmr_data-class}} containing
#'   the SNP-level summary statistics.
#' @param random_start A length-one logical. If \code{TRUE} (default),
#'   starting values are drawn randomly from diffuse distributions. If
#'   \code{FALSE}, the values in \code{start_values} are validated and
#'   returned.
#' @param K_start A length-one positive integer giving the initial number of
#'   mixture clusters. Use \code{K_start = 1} for the non-mixture models
#'   (\code{\link{bayesmr}()}, \code{\link{bayesmr_het}()}).
#' @param start_values A named list of user-supplied starting values, used
#'   when \code{random_start = FALSE}. See \code{\link{bayesmr_startvalues}()}
#'   for the expected format.
#'
#' @return A named list of starting values. For \code{K_start = 1}: contains
#'   \code{gamma}, \code{beta}, \code{psi}, \code{tau}. For
#'   \code{K_start > 1}: also contains \code{K}, \code{xi}, \code{alpha}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso \code{\link{bayesmr}()}, \code{\link{bayesmr_startvalues}()},
#'   \code{\link{bayesmr_data-class}}.
#'
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical
#'   Modeling for Two-Sample Summary-Data Mendelian Randomization under
#'   Heterogeneity", working paper.
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' # Random start for the no-clustering model
#' sv1 <- bayesmr_init(dat, random_start = TRUE, K_start = 1)
#' # Random start for the mixture model with 3 initial clusters
#' sv3 <- bayesmr_init(dat, random_start = TRUE, K_start = 3)
#' }
#'
#' @export
bayesmr_init <- function(data, random_start = TRUE, K_start = 1, start_values) {
  if (random_start) {
    # initialize psi and tau
    psi_tau <- rhalft(2, alpha = .01, nu = 3)
    psi <- psi_tau[1]
    tau <- psi_tau[2]

    if (K_start == 1) {
      # initialize gamma and beta
      gamma_beta <- rnorm(2, mean = 0, sd = .5)
      gamma <- gamma_beta[1]
      beta <- gamma_beta[2]

      return(list(gamma = gamma, beta = beta, psi = psi, tau = tau))
    }
    else if (K_start > 1) {
      # initialize gamma
      gamma <- rnorm(1, mean = 0, sd = .5)

      # initialize xi (cluster indicators)
      xi <- sample(x = 1:K_start, size = data@n, replace = TRUE)
      while (length(table(xi)) != K_start) {  # make sure the K_start clusters are all non empty
        xi <- sample(x = 1:K_start, size = data@n, replace = TRUE)
      }

      # initialize beta unique values
      beta <- rnorm(K_start, 0, sd = .5)      # beta_star
    }
      
    # initialize alpha (concentration)
    alpha <- 1

    return(list(gamma = gamma, beta = beta, psi = psi, tau = tau, K = K_start, xi = xi, alpha = alpha))
  }
  else {
    if (check_startvalues(start_values)) {
      return(start_values)
    }
    else {
      stop("the starting values list is not correct; see the documentation for more details.")
    }
  }
}

#' Construct a Starting-Values List for BayesMR
#'
#' \code{bayesmr_startvalues()} assembles a named list of user-supplied
#' starting values to pass to \code{\link{bayesmr_init}()} when
#' \code{random_start = FALSE}.
#'
#' \code{startvalues_bayesmr()} is an alias for \code{bayesmr_startvalues()}.
#'
#' \code{check_startvalues()} validates the list before model fitting.
#'
#' @param gamma A length-one numeric value for the initial causal effect
#'   \eqn{\gamma}.
#' @param beta A numeric vector of initial causal effect(s) \eqn{\beta} (or
#'   cluster-specific effects \eqn{\beta^*_1, \ldots, \beta^*_K} for mixture
#'   models). Length must equal \code{K} (supplied via \code{...}).
#' @param ... Additional named starting values. For mixture models, supply
#'   \code{K} (integer, number of clusters), \code{xi} (integer vector of
#'   length \eqn{n} with cluster labels in \eqn{1,\ldots,K}), and
#'   \code{alpha} (positive scalar, DP concentration). For heterogeneity
#'   models, also supply \code{psi} and \code{tau} (positive scalars).
#' @param startvalues A named list of starting values (used by
#'   \code{check_startvalues()}).
#'
#' @return \code{bayesmr_startvalues()} returns a named list.
#'   \code{check_startvalues()} returns a length-one logical: \code{TRUE} if
#'   the list is valid, \code{FALSE} otherwise.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso \code{\link{bayesmr_init}()}, \code{\link{bayesmr}()}.
#'
#' @examples
#' # Starting values for the no-clustering model
#' sv <- bayesmr_startvalues(gamma = 0.1, beta = 0.2,
#'                           K = 1, psi = 0.05, tau = 0.05)
#' check_startvalues(sv)
#'
#' @export
bayesmr_startvalues <- function(gamma, beta, ...) {
  c(list(gamma = gamma, beta = beta), list(...))
}

#' @rdname bayesmr_startvalues
#' @export
startvalues_bayesmr <- bayesmr_startvalues

#' @rdname bayesmr_startvalues
#' @export
check_startvalues <- function(startvalues) {
  startvalues_ok <- TRUE

  # check startvalues list
  if (!is.list(startvalues)) {
    startvalues_ok <- FALSE
    return(startvalues_ok)
  }

  sv_names <- names(startvalues)

  # check gamma startvalues
  if (length(startvalues[["gamma"]]) != 1) {
    startvalues_ok <- FALSE
    return(startvalues_ok)
  }
  # check beta startvalues
  if (is.null(startvalues[["beta"]])) {
    startvalues_ok <- FALSE
    return(startvalues_ok)
  }

  # check K startvalues
  if (!is.null(startvalues[["K"]])) {
    if (length(startvalues[["K"]]) != 1) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }
    if (startvalues[["K"]] < 1) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }
    if (startvalues[["K"]] != trunc(startvalues[["K"]])) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }

    if (length(startvalues[["beta"]]) != startvalues[["K"]]) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }

    if (startvalues[["K"]] > 1) {
      # check xi startvalues
      if (is.null(startvalues[["xi"]])) {
        startvalues_ok <- FALSE
        return(startvalues_ok)
      }
      if (any(startvalues[["xi"]] != trunc(startvalues[["xi"]]))) {
        startvalues_ok <- FALSE
        return(startvalues_ok)
      }
      if (any(startvalues[["xi"]] < 1)) {
        startvalues_ok <- FALSE
        return(startvalues_ok)
      }

      # check alpha startvalues
      if (is.null(startvalues[["alpha"]])) {
        startvalues_ok <- FALSE
        return(startvalues_ok)
      }
      if (length(startvalues[["alpha"]]) != 1) {
        startvalues_ok <- FALSE
        return(startvalues_ok)
      }
      if (startvalues[["alpha"]] <= 0) {
        startvalues_ok <- FALSE
        return(startvalues_ok)
      }
    }
  }
  else {
    startvalues_ok <- FALSE
    return(startvalues_ok)
  }

  # check psi startvalues
  if (!is.null(startvalues[["psi"]])) {
    if (length(startvalues[["psi"]]) != 1) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }
    if (startvalues[["psi"]] < 0) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }
  }

  # check tau startvalues
  if (!is.null(startvalues[["tau"]])) {
    if (length(startvalues[["tau"]]) != 1) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }
    if (startvalues[["tau"]] < 0) {
      startvalues_ok <- FALSE
      return(startvalues_ok)
    }
  }

  return(startvalues_ok)
}
