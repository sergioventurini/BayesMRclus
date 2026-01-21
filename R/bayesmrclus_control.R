#' Auxiliary Function for Controlling BayesM Model Fitting
#' 
#' @description{
#' \code{bayesmr_control()} is an auxiliary function as user interface for
#'   \code{bayesmr()} fitting. Typically only used when calling the \code{bayesmr()}
#'   function. It is used to set parameters that affect the sampling but do
#'   not affect the posterior distribution.
#' 
#' \code{control_bayesmr()} is an alias for \code{bayesmr_control()}.
#' 
#' \code{check_control()} is an auxiliary function that verifies the
#'   correctness of the controls provided before a BayesM is fitted with
#'   \code{\link{bayesmr}()}.
#' }
#'
#' @param nsim A length-one numeric vector for the number of draws to be taken
#'   from the posterior distribution.
#' @param burnin A length-one numeric vector for the number of initial MCMC
#'   iterations (usually to be discarded).
#' @param thin A length-one numeric vector for the number of iterations between
#'   consecutive draws.
#' @param nchains A length-one numeric vector for the number of parallel chains to run.
#' @param threads A length-one numeric vector for the number of chains to run.
#'   If greater than 1, package \pkg{\link{parallel}} is used to take advantage of any
#'   multiprocessing or distributed computing capabilities that may be available.
#' @param seed An integer scalar. If supplied, provides the random number seed.
#' @param parallel A length-one character vector indicating the type of parallel
#'   operation to be used (if any). Possible values are \code{multicore}
#'   (which works only on Unix/mcOS), \code{snow} and \code{no} (i.e. serial
#'   instead of parallel computing).
#' @param beta.prop A length-one numeric vector providing the standard deviation of the
#'   proposal distribution for the jump in the individual latent space
#'   position.
#' @param beta.m A length-one integer vector providing the number of beta values to draw.
#' @param random.start A length-one logical vector. If \code{TRUE} the starting
#'   values are drawn randomly, otherwise a user-defined starting partition must
#'   be provided through the \code{partition} argument.
#' @param partition A length-one numeric vector providing the user-defined
#'   starting partition.
#' @param procrustes A length-one logical vector. If \code{TRUE} the simulated
#'   MCMC chains are post-processed through a Procrustes transformation.
#' @param relabel A length-one logical vector. If \code{TRUE} the simulated
#'   MCMC chains are relabelled to address the label-switching problem.
#' @param store.burnin A logical scalar. If \code{TRUE}, the samples from the
#'   burnin are also stored and returned.
#' @param verbose A logical scalar. If \code{TRUE}, causes information to be
#'   printed out about the progress of the fitting.
#' @param control A list of control options.
#' @return A named list with the control options as components.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bayesmr}()}
#' @keywords model based clustering
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' # Shorter run than default.
#' sim.fit <- bayesmr(simdiss,
#'   control = bayesmr_control(burnin = 1000, nsim = 2000, thin = 5, verbose = TRUE))
#' }
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
                            random.start = TRUE,
                            partition = NULL,
                            procrustes = TRUE,
                            relabel = TRUE,
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
  if (!is.logical(control[["random.start"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!control[["random.start"]]) {
    if (!is.numeric(control[["partition"]])) {
      control_ok <- FALSE
      return(control_ok)
    }
  }
  if (!is.logical(control[["store.burnin"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["verbose"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["procrustes"]])) {
    control_ok <- FALSE
    return(control_ok)
  }
  if (!is.logical(control[["relabel"]])) {
    control_ok <- FALSE
    return(control_ok)
  }

  return(control_ok)
}
