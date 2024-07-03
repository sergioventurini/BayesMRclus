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
#' @param eta A named list containing the hyperparameters for the prior
#'   distribution of the \eqn{\eta_1,\ldots,\eta_G} parameters. It should
#'   contain two numeric vectors, namely \code{a} and \code{b}.
#' @param sigma2 A named list containing the hyperparameters for the prior
#'   distributions of the \eqn{\sigma^2_1,\ldots,\sigma^2_G} parameters. It
#'   should contain two numeric scalars, namely \code{a} and \code{b}.
#' @param lambda A list containing the hyperparameters for the prior
#'   distribution of the \eqn{\lambda_1,\ldots,\lambda_G} parameters. It should
#'   contain a single numeric vector.
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
#'   prior = bayesmr_prior(sigma2 = list(a = 1, b = 4)))
#' }
#'
#' @export
bayesmr_prior <- function(eta = list(a = rep(1.5, .bayesmrEnv$current_G), b = rep(.5, .bayesmrEnv$current_G)),
                       sigma2 = list(a = 1e-1, b = 1e-1),
                       lambda = rep(1, .bayesmrEnv$current_G)){
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

  # check eta prior
  if (!is.list(prior[["eta"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["eta"]][["a"]]) != .bayesmrEnv$current_G) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["eta"]][["b"]]) != .bayesmrEnv$current_G) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["eta"]][["a"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["eta"]][["b"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check sigma2 prior
  if (!is.list(prior[["sigma2"]])) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (length(prior[["sigma2"]]) != 2) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["sigma2"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  # check lambda prior
  if (length(prior[["lambda"]]) != .bayesmrEnv$current_G) {
    prior_ok <- FALSE
    return(prior_ok)
  }
  if (any(prior[["lambda"]] < 0)) {
    prior_ok <- FALSE
    return(prior_ok)
  }

  return(prior_ok)
}
