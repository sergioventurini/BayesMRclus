#' Function to compute the starting values before fitting a BayesMR models.
#'
#' \code{bayesmr_init()} is the main function that estimates a BayesMR model.
#'
#' @param data A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param random_start A length-one logical vector. If \code{TRUE} the starting
#'   values are drawn randomly, otherwise.
#' @param K_start A length-one numeric vector providing the initial
#'   value for the number of clusters.
#' @return A named \code{list} with the following items:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates starting values}
#'     \item{\code{x}: }{numeric vector of initial cluster memberships}
#'     \item{\code{ng}: }{numeric vector of initial cluster sizes}
#'     \item{\code{alpha}: }{numeric vector of alpha starting values}
#'     \item{\code{eta}: }{numeric vector of eta starting values}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 starting values}
#'     \item{\code{lambda}: }{numeric vector of lambda starting values}
#'   }
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bayesmr}()} for fitting a BayesMR model.
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#' @examples
#' data(simdiss, package = "bayesmr")
#' bayesmr_init(simdiss@diss, random_start = TRUE)
#' @export
bayesmr_init <- function(data, random_start = TRUE, K_start = 1, start_values) {
  n <- nrow(data)
  
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

#' Auxiliary Function for Setting BayesMR Model Starting Values
#' 
#' @description{
#' \code{bayesmr_startvalues()} is an auxiliary function as user interface for
#'   \code{bayesmr()} fitting. Typically only used when calling the \code{bayesmr()}
#'   function. It is used to set prior hyperparameters.
#' 
#' \code{prior_bayesmr()} is an alias for \code{bayesmr_startvalues()}.
#' 
#' \code{check_prior()} is an auxiliary function that verifies the
#'   correctness of the prior hyperparameters provided before a BayesMR is fitted
#'   with \code{\link{bayesmr}()}.
#' 
#' \code{update_prior()} is an auxiliary function to modify a set of prior
#'   choices using a new value of \emph{p} and \emph{G}. It is intended for
#'   internal use mainly in the \code{\link{bayesmr_ic}()} function.
#' }
#'
#' @param gammaj A named list containing the hyperparameters for the prior
#'   distribution of the \eqn{\eta_1,\ldots,\eta_G} parameters. It should
#'   contain two numeric vectors, namely \code{a} and \code{b}.
#' @param Gammaj A named list containing the hyperparameters for the prior
#'   distribution of the \eqn{\eta_1,\ldots,\eta_G} parameters. It should
#'   contain two numeric vectors, namely \code{a} and \code{b}.
#' @param gamma A named list containing the hyperparameters for the prior
#'   distribution of the \eqn{\eta_1,\ldots,\eta_G} parameters. It should
#'   contain two numeric vectors, namely \code{a} and \code{b}.
#' @param beta A named list containing the hyperparameters for the prior
#'   distributions of the \eqn{\sigma^2_1,\ldots,\sigma^2_G} parameters. It
#'   should contain two numeric scalars, namely \code{a} and \code{b}.
#' @param prior A named list of prior hyperparameters.
#' @return A list with the prior hyperparameters as components.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bayesmr}()}
#' @keywords model based clustering
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' # Shorter run than default.
#' sim.fit <- bayesmr(simdiss,
#'   control = bayesmr_control(burnin = 1000, nsim = 2000, thin = 1, verbose = TRUE),
#'   prior = bayesmr_startvalues(gamma = list(mean = 0, var = 1)))
#' }
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
  if (startvalues[["gamma"]] < 0) {
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
      if (length(startvalues[["xi"]]) != startvalues[["K"]]) {
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
