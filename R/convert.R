#' Conversion of an \code{bayesmr_fit} object to an object of class \code{mcmc}.
#' 
#' \code{bayesmr_fit_to_mcmc} converts an object of class \code{bayesmr_fit}
#'   to an object with class \code{mcmc}.
#' 
#' @param res An object of type \code{bayesmr_fit}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc}.
#' @seealso
#'   \code{\link{bayesmr}()} for for fitting a DMBC model;
#'   \code{\link{bayesmr_fit-class}};
#'   \code{\link[coda]{mcmc}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], verbose = TRUE)
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#' sim.mcmc <- bayesmr_fit_to_mcmc(sim.bayesmr@results[[1]], TRUE)
#' plot(sim.mcmc)
#' }
#' @export
bayesmr_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@dim[["n"]]
  p <- res@dim[["p"]]
  G <- res@dim[["G"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      warning("burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  theta <- matrix(NA, nrow = length(tokeep), ncol = 2 + 3)
  theta_nm <- character(2 + 3)
  theta[, 1] <- res@gamma.chain[tokeep]
  theta_nm[1] <- "gamma"
  theta[, 2] <- res@beta.chain[tokeep]
  theta_nm[2] <- "beta"
  theta[, 2 + 1] <- res@dens$loglik[tokeep]
  theta_nm[2 + 1] <- "loglik"
  theta[, 2 + 2] <- res@dens$logprior[tokeep]
  theta_nm[2 + 2] <- "logprior"
  theta[, 2 + 3] <- res@dens$logpost[tokeep]
  theta_nm[2 + 3] <- "logpost"
  colnames(theta) <- theta_nm

  if (store.burnin) {
    if (include.burnin) {
      out <- coda::mcmc(theta, start = 1, end = totiter, thin = thin)
    } else {
      out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
  }

  return(out)
}

#' Conversion of an \code{bayesmr_fit_list} object to a \code{list}.
#' 
#' \code{bayesmr_fit_list_to_list} converts an object of class
#'   \code{bayesmr_fit_list} to a list of arrays including all the parameter.
#'   chains. It is intended for internal use mainly.
#' 
#' @param res An object of type \code{bayesmr_fit_list}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc.list}.
#' @seealso
#'   \code{\link{bayesmr}()} for for fitting a DMBC model;
#'   \code{\link{bayesmr_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], nchains = 2, verbose = TRUE)
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#' sim.list <- bayesmr_fit_list_to_list(sim.bayesmr, TRUE)
#' 
#' library(bayesplot)
#' mcmc_trace(sim.list, regex_pars = "lambda")
#' }
#' @export
bayesmr_fit_list_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@results[[1]]@dim[["n"]]
  p <- res@results[[1]]@dim[["p"]]
  G <- res@results[[1]]@dim[["G"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      warning("burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  out <- list()
  for (c in 1:nchains) {
    theta <- matrix(NA, nrow = length(tokeep), ncol = (2 + 3))
    theta_nm <- character(2 + 3)
    theta[, 1] <- res@results[[c]]@gamma.chain[tokeep]
    theta_nm[1] <- "gamma"
    theta[, 2] <- res@results[[c]]@beta.chain[tokeep]
    theta_nm[2] <- "beta"

    theta[, 2 + 1] <- res@results[[c]]@dens$loglik[tokeep]
    theta_nm[2 + 1] <- "loglik"
    theta[, 2 + 2] <- res@results[[c]]@dens$logprior[tokeep]
    theta_nm[2 + 2] <- "logprior"
    theta[, 2 + 3] <- res@results[[c]]@dens$logpost[tokeep]
    theta_nm[2 + 3] <- "logpost"
    colnames(theta) <- theta_nm
    out[[c]] <- theta
  }

  return(out)
}

#' Conversion of an \code{bayesmr_fit_list} object to an object of class
#'   \code{mcmc.list}.
#' 
#' \code{bayesmr_fit_list_to_mcmc.list} converts an object of class
#'   \code{bayesmr_fit_list} to an object with class \code{mcmc.list}.
#' 
#' @param res An object of type \code{bayesmr_fit_list}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc.list}.
#' @seealso
#'   \code{\link{bayesmr}()} for for fitting a DMBC model;
#'   \code{\link{bayesmr_fit_list-class}};
#'   \code{\link[coda]{mcmc.list}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], nchains = 2, verbose = TRUE)
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#' sim.mcmc <- bayesmr_fit_list_to_mcmc.list(sim.bayesmr, TRUE)
#' plot(sim.mcmc)
#' }
#' @export
bayesmr_fit_list_to_mcmc.list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  out <- bayesmr_fit_list_to_list(res, include.burnin = include.burnin, verbose = verbose)
  if (store.burnin) {
    if (include.burnin) {
      out <- lapply(out, coda::mcmc, start = 1, end = totiter, thin = thin)
    } else {
      out <- lapply(out, coda::mcmc, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- lapply(out, coda::mcmc, start = (burnin + 1), end = totiter, thin = thin)
  }
  out <- coda::mcmc.list(out)

  return(out)
}
































#' Conversion of an \code{bayesmr_full_fit} object to an object of class \code{mcmc}.
#' 
#' \code{bayesmr_full_fit_to_mcmc} converts an object of class \code{bayesmr_full_fit}
#'   to an object with class \code{mcmc}.
#' 
#' @param res An object of type \code{bayesmr_full_fit}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc}.
#' @seealso
#'   \code{\link{bayesmr}()} for for fitting a DMBC model;
#'   \code{\link{bayesmr_full_fit-class}};
#'   \code{\link[coda]{mcmc}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], verbose = TRUE)
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#' sim.mcmc <- bayesmr_full_fit_to_mcmc(sim.bayesmr@results[[1]], TRUE)
#' plot(sim.mcmc)
#' }
#' @export
bayesmr_full_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@dim[["n"]]
  p <- res@dim[["p"]]
  G <- res@dim[["G"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      warning("burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  theta <- matrix(NA, nrow = length(tokeep), ncol = 4 + 3)
  theta_nm <- character(4 + 3)
  theta[, 1] <- res@gamma.chain[tokeep]
  theta_nm[1] <- "gamma"
  theta[, 2] <- res@beta.chain[tokeep]
  theta_nm[2] <- "beta"
  theta[, 3] <- res@psi.chain[tokeep]
  theta_nm[3] <- "psi"
  theta[, 4] <- res@tau.chain[tokeep]
  theta_nm[4] <- "tau"
  theta[, 4 + 1] <- res@dens$loglik[tokeep]
  theta_nm[4 + 1] <- "loglik"
  theta[, 4 + 2] <- res@dens$logprior[tokeep]
  theta_nm[4 + 2] <- "logprior"
  theta[, 4 + 3] <- res@dens$logpost[tokeep]
  theta_nm[4 + 3] <- "logpost"
  colnames(theta) <- theta_nm

  if (store.burnin) {
    if (include.burnin) {
      out <- coda::mcmc(theta, start = 1, end = totiter, thin = thin)
    } else {
      out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
  }

  return(out)
}

#' Conversion of an \code{bayesmr_full_fit_list} object to a \code{list}.
#' 
#' \code{bayesmr_full_fit_list_to_list} converts an object of class
#'   \code{bayesmr_full_fit_list} to a list of arrays including all the parameter.
#'   chains. It is intended for internal use mainly.
#' 
#' @param res An object of type \code{bayesmr_full_fit_list}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc.list}.
#' @seealso
#'   \code{\link{bayesmr}()} for for fitting a DMBC model;
#'   \code{\link{bayesmr_full_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], nchains = 2, verbose = TRUE)
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#' sim.list <- bayesmr_full_fit_list_to_list(sim.bayesmr, TRUE)
#' 
#' library(bayesplot)
#' mcmc_trace(sim.list, regex_pars = "lambda")
#' }
#' @export
bayesmr_full_fit_list_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@results[[1]]@dim[["n"]]
  p <- res@results[[1]]@dim[["p"]]
  G <- res@results[[1]]@dim[["G"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      warning("burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  out <- list()
  for (c in 1:nchains) {
    theta <- matrix(NA, nrow = length(tokeep), ncol = (4 + 3))
    theta_nm <- character(4 + 3)
    theta[, 1] <- res@results[[c]]@gamma.chain[tokeep]
    theta_nm[1] <- "gamma"
    theta[, 2] <- res@results[[c]]@beta.chain[tokeep]
    theta_nm[2] <- "beta"
    theta[, 3] <- res@results[[c]]@psi.chain[tokeep]
    theta_nm[3] <- "psi"
    theta[, 4] <- res@results[[c]]@tau.chain[tokeep]
    theta_nm[4] <- "tau"

    theta[, 4 + 1] <- res@results[[c]]@dens$loglik[tokeep]
    theta_nm[4 + 1] <- "loglik"
    theta[, 4 + 2] <- res@results[[c]]@dens$logprior[tokeep]
    theta_nm[4 + 2] <- "logprior"
    theta[, 4 + 3] <- res@results[[c]]@dens$logpost[tokeep]
    theta_nm[4 + 3] <- "logpost"
    colnames(theta) <- theta_nm
    out[[c]] <- theta
  }

  return(out)
}

#' Conversion of an \code{bayesmr_full_fit_list} object to an object of class
#'   \code{mcmc.list}.
#' 
#' \code{bayesmr_full_fit_list_to_mcmc.list} converts an object of class
#'   \code{bayesmr_full_fit_list} to an object with class \code{mcmc.list}.
#' 
#' @param res An object of type \code{bayesmr_full_fit_list}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc.list}.
#' @seealso
#'   \code{\link{bayesmr}()} for for fitting a DMBC model;
#'   \code{\link{bayesmr_full_fit_list-class}};
#'   \code{\link[coda]{mcmc.list}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], nchains = 2, verbose = TRUE)
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#' sim.mcmc <- bayesmr_full_fit_list_to_mcmc.list(sim.bayesmr, TRUE)
#' plot(sim.mcmc)
#' }
#' @export
bayesmr_full_fit_list_to_mcmc.list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  out <- bayesmr_full_fit_list_to_list(res, include.burnin = include.burnin, verbose = verbose)
  if (store.burnin) {
    if (include.burnin) {
      out <- lapply(out, coda::mcmc, start = 1, end = totiter, thin = thin)
    } else {
      out <- lapply(out, coda::mcmc, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- lapply(out, coda::mcmc, start = (burnin + 1), end = totiter, thin = thin)
  }
  out <- coda::mcmc.list(out)

  return(out)
}
