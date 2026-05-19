#' Auxiliary Function for Controlling BayesMR Model Fitting
#'
#' \code{bayesmr_control()} creates the list of MCMC control parameters used
#' by \code{\link{bayesmr}()} and its variants. It sets parameters that affect
#' the sampling algorithm but not the posterior distribution.
#'
#' \code{control_bayesmr()} is an alias for \code{bayesmr_control()}.
#'
#' \code{check_control()} validates a control list before model fitting.
#'
#' @param nsim A length-one positive integer giving the number of post-burn-in
#'   MCMC draws to take from the posterior distribution.
#' @param burnin A length-one non-negative integer giving the number of
#'   initial MCMC iterations to discard as burn-in.
#' @param thin A length-one positive integer giving the thinning interval;
#'   only every \code{thin}-th draw is retained.
#' @param nchains A length-one positive integer giving the number of
#'   independent MCMC chains to run.
#' @param threads A length-one positive integer giving the number of parallel
#'   workers to use. If greater than 1 and \code{parallel != "no"},
#'   \pkg{parallel} is used to distribute chains across workers.
#' @param seed An optional integer scalar providing the random number seed.
#'   If \code{NULL} (default), the current RNG state is used.
#' @param parallel A length-one character vector indicating the parallelism
#'   backend: \code{"no"} (serial, default), \code{"multicore"} (Unix/macOS
#'   only, via \code{\link[parallel]{mclapply}}), or \code{"snow"} (via
#'   \code{\link[parallel]{parLapply}}).
#' @param beta.prop A length-one positive numeric giving the standard deviation
#'   of the random-walk Metropolis-Hastings proposal for \eqn{\beta}.
#' @param beta.m A length-one positive integer giving the number of \eqn{\beta}
#'   values drawn per mixture cluster at each iteration.
#' @param psi.prop A length-one positive numeric giving the standard deviation
#'   of the random-walk Metropolis-Hastings proposal for \eqn{\psi}.
#' @param tau.prop A length-one positive numeric giving the standard deviation
#'   of the random-walk Metropolis-Hastings proposal for \eqn{\tau}.
#' @param random_start A logical scalar. If \code{TRUE} (default), starting
#'   values are drawn randomly via \code{\link{bayesmr_init}()}; otherwise,
#'   user-supplied starting values are used.
#' @param K_start A length-one positive integer giving the initial number of
#'   mixture clusters.
#' @param store.burnin A logical scalar. If \code{TRUE} (default), burn-in
#'   iterations are stored and returned (useful for diagnosing convergence).
#' @param verbose A logical scalar. If \code{TRUE}, prints progress information
#'   during fitting.
#' @param control A named list of control options (used by \code{check_control()}).
#'
#' @return \code{bayesmr_control()} (and its alias) returns a named list with
#'   all control parameters as components.
#'   \code{check_control()} returns a length-one logical: \code{TRUE} if the
#'   list is valid, \code{FALSE} otherwise.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso \code{\link{bayesmr}()}, \code{\link{bayesmr_prior}()}.
#'
#' @examples
#' # Default control settings
#' ctrl <- bayesmr_control()
#' str(ctrl)
#'
#' # Custom settings for a short run
#' ctrl2 <- bayesmr_control(nsim = 1000, burnin = 500, nchains = 2,
#'                          threads = 2, parallel = "multicore")
#'
#' @export
bayesmr_control <- function(nsim = 5000,
                            burnin = 10000,
                            thin = 1,
                            nchains = 1,
                            threads = 1,
                            seed = NULL,
                            parallel = "no",
                            beta.prop = .5,
                            beta.m = 2,
                            psi.prop = .1,
                            tau.prop = .1,
                            random_start = TRUE,
                            K_start = 2,
                            store.burnin = TRUE,
                            verbose = FALSE){
  control <- list()
  for (arg in names(formals(sys.function())))
    control[[arg]] <- get(arg)
  control
}

#' @rdname bayesmr_control
#' @export
control_bayesmr <- bayesmr_control

#' @rdname bayesmr_control
#' @export
check_control <- function(control) {
  control_ok <- TRUE

  if (!is.list(control)) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["nsim"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["burnin"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["thin"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["nchains"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["threads"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.null(control[["seed"]])) {
    if (control[["seed"]] < 1) {
      control_ok <- FALSE
      return(control_ok)
    }
  }
  if (!(control[["parallel"]] %in% c("no", "snow", "multicore"))) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["beta.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["beta.m"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["psi.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["tau.prop"]] < 0) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["random_start"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.numeric(control[["K_start"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (control[["K_start"]] < 1) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["store.burnin"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["verbose"]])) {
    control_ok <- FALSE
    return(control_ok)
  }

  return(control_ok)
}
