#' Auxiliary Function for Setting BayesMR Prior Hyperparameters
#'
#' \code{bayesmr_prior()} creates the named list of prior hyperparameters used
#' by \code{\link{bayesmr}()} and its variants.
#'
#' \code{prior_bayesmr()} is an alias for \code{bayesmr_prior()}.
#'
#' \code{check_prior()} validates a prior list before model fitting.
#'
#' @param gammaj A named list with element \code{psi2}: the fixed variance
#'   component \eqn{\psi^2} added to each SNP-exposure variance
#'   \eqn{\sigma^2_{X,j}} to form \eqn{\psi^2_j = \sigma^2_{X,j} + \psi^2}.
#'   Used only in the fixed-heterogeneity models (\code{\link{bayesmr}()},
#'   \code{\link{bayesmr_mix}()}).
#' @param Gammaj A named list with element \code{tau2}: the fixed variance
#'   component \eqn{\tau^2} added to each SNP-outcome variance
#'   \eqn{\sigma^2_{Y,j}} to form \eqn{\tau^2_j = \sigma^2_{Y,j} + \tau^2}.
#'   Used only in the fixed-heterogeneity models.
#' @param psi A named list with elements \code{alpha} (scale) and \code{nu}
#'   (degrees of freedom) of the half-\eqn{t} prior on \eqn{\psi}
#'   (exposure heterogeneity scale). Used in \code{\link{bayesmr_het}()} and
#'   \code{\link{bayesmr_mix_het}()}.
#' @param tau A named list with elements \code{alpha} (scale) and \code{nu}
#'   (degrees of freedom) of the half-\eqn{t} prior on \eqn{\tau}
#'   (outcome heterogeneity scale). Used in \code{\link{bayesmr_het}()} and
#'   \code{\link{bayesmr_mix_het}()}.
#' @param gamma A named list with elements \code{mean} and \code{var} for the
#'   Normal prior on the global causal effect \eqn{\gamma}.
#' @param beta A named list with elements \code{mean} and \code{var} for the
#'   Normal prior on \eqn{\beta} (or on each cluster-specific causal effect
#'   in mixture models).
#' @param alpha A named list with elements \code{a} and \code{b} for the
#'   Beta(\eqn{a}, \eqn{b}) prior on the Dirichlet Process concentration
#'   parameter \eqn{\alpha}. Used in \code{\link{bayesmr_mix}()} and
#'   \code{\link{bayesmr_mix_het}()}.
#' @param prior A named list of prior hyperparameters (used by
#'   \code{check_prior()}).
#'
#' @return \code{bayesmr_prior()} (and its alias) returns a named list with
#'   all hyperparameters as components. \code{check_prior()} returns a
#'   length-one logical: \code{TRUE} if the list is valid, \code{FALSE}
#'   otherwise.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso \code{\link{bayesmr}()}, \code{\link{bayesmr_control}()},
#'   \code{\link{dhalft}()}.
#'
#' @examples
#' # Default prior hyperparameters
#' pr <- bayesmr_prior()
#' str(pr)
#'
#' # Tighter normal priors on gamma and beta
#' pr2 <- bayesmr_prior(gamma = list(mean = 0, var = 0.25),
#'                      beta  = list(mean = 0, var = 0.25))
#'
#' @export
bayesmr_prior <- function(gammaj = list(psi2 = 1), Gammaj = list(tau2 = 1),
                          psi = list(alpha = 0.1, nu = 3),
                          tau = list(alpha = 0.1, nu = 3),
                          gamma = list(mean = 0, var = 1),
                          beta = list(mean = 0, var = 1),
                          alpha = list(a = 2, b = 2)){
  prior <- list()
  for (arg in names(formals(sys.function())))
    prior[[arg]] <- get(arg)
  prior
}

#' @rdname bayesmr_prior
#' @export
prior_bayesmr <- bayesmr_prior

#' @rdname bayesmr_prior
#' @export
check_prior <- function(prior) {
  prior_ok <- TRUE

  # check prior list
  if (!is.list(prior)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check gammaj prior
  if (!is.list(prior[["gammaj"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["gammaj"]]) != 1) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["gammaj"]][["psi2"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check Gammaj prior
  if (!is.list(prior[["Gammaj"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["Gammaj"]]) != 1) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["Gammaj"]][["tau2"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check psi prior
  if (!is.list(prior[["psi"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["psi"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["psi"]][["alpha"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["psi"]][["nu"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check tau prior
  if (!is.list(prior[["tau"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["tau"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["tau"]][["alpha"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["tau"]][["nu"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check gamma prior
  if (!is.list(prior[["gamma"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["gamma"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["gamma"]]["var"] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check beta prior
  if (!is.list(prior[["beta"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["beta"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["beta"]]["var"] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check alpha prior
  if (!is.list(prior[["alpha"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["alpha"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["alpha"]] <= 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  return(prior_ok)
}
