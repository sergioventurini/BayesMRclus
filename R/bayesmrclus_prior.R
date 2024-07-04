#' Auxiliary Function for Setting BayesMR Model Priors
#' 
#' @description{
#' \code{bayesmr_prior()} is an auxiliary function as user interface for
#'   \code{bayesmr()} fitting. Typically only used when calling the \code{bayesmr()}
#'   function. It is used to set prior hyperparameters.
#' 
#' \code{prior_bayesmr()} is an alias for \code{bayesmr_prior()}.
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
#'   prior = bayesmr_prior(gamma = list(mean = 0, var = 1)))
#' }
#'
#' @export
bayesmr_prior <- function(gammaj = list(psi2 = 1), Gammaj = list(tau2 = 1),
                          gamma = list(mean = 0, var = 1),
                          beta = list(mean = 0, var = 1)){
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
  if (any(prior[["gammaj"]][["psi2"]] <= 0)) {
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
  if (any(prior[["Gammaj"]][["tau2"]] <= 0)) {
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

  return(prior_ok)
}
