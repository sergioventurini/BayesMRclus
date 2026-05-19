#' Conversion of a \code{bayesmr_fit} Object to \code{mcmc}
#'
#' Converts a single-chain \code{\link{bayesmr_fit-class}} object to a
#' \code{\link[coda]{mcmc}} object containing the \eqn{\gamma}, \eqn{\beta},
#' and log-density chains.
#'
#' @param res An object of class \code{bayesmr_fit}.
#' @param include.burnin Logical; if \code{TRUE} the burn-in iterations
#'   (when stored) are included. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} prints warnings during conversion.
#'   Default is \code{TRUE}.
#'
#' @return A \code{\link[coda]{mcmc}} object.
#'
#' @seealso \code{\link{bayesmr}()}, \code{\link{bayesmr_fit-class}},
#'   \code{\link[coda]{mcmc}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' res <- bayesmr(dat, control = bayesmr_control(nsim = 500, burnin = 500))
#' mc <- bayesmr_fit_to_mcmc(res@results[[1]])
#' plot(mc)
#' }
#'
#' @export
bayesmr_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@dim[["n"]]

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

#' Conversion of a \code{bayesmr_fit_list} Object to a \code{list}
#'
#' Converts a multi-chain \code{\link{bayesmr_fit_list-class}} object to a
#' plain list of matrices (one per chain) suitable for use with
#' \pkg{bayesplot} functions.
#'
#' @param res An object of class \code{bayesmr_fit_list}.
#' @param include.burnin Logical; if \code{TRUE} the burn-in iterations
#'   (when stored) are included. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} prints warnings during conversion.
#'   Default is \code{TRUE}.
#'
#' @return A list of numeric matrices, one per chain, with columns
#'   \code{gamma}, \code{beta}, \code{loglik}, \code{logprior},
#'   \code{logpost}.
#'
#' @seealso \code{\link{bayesmr}()}, \code{\link{bayesmr_fit_list-class}},
#'   \code{\link{bayesmr_fit_list_to_mcmc.list}()}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' res <- bayesmr(dat, control = bayesmr_control(nsim = 500, burnin = 500,
#'                                               nchains = 2))
#' lst <- bayesmr_fit_list_to_list(res)
#' bayesplot::mcmc_trace(lst, pars = "gamma")
#' }
#'
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

#' Conversion of a \code{bayesmr_fit_list} Object to \code{mcmc.list}
#'
#' Converts a multi-chain \code{\link{bayesmr_fit_list-class}} object to a
#' \code{\link[coda]{mcmc.list}} for convergence diagnostics.
#'
#' @param res An object of class \code{bayesmr_fit_list}.
#' @param include.burnin Logical; if \code{TRUE} the burn-in iterations
#'   (when stored) are included. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} prints warnings during conversion.
#'   Default is \code{TRUE}.
#'
#' @return A \code{\link[coda]{mcmc.list}} object.
#'
#' @seealso \code{\link{bayesmr}()}, \code{\link{bayesmr_fit_list-class}},
#'   \code{\link[coda]{mcmc.list}}, \code{\link{bayesmr_fit_list_to_list}()}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' res <- bayesmr(dat, control = bayesmr_control(nsim = 500, burnin = 500,
#'                                               nchains = 2))
#' mc <- bayesmr_fit_list_to_mcmc.list(res)
#' coda::gelman.diag(mc)
#' }
#'
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

#' Conversion of a \code{bayesmr_het_fit} Object to \code{mcmc}
#'
#' Converts a single-chain \code{\link{bayesmr_het_fit-class}} object to a
#' \code{\link[coda]{mcmc}} object containing the \eqn{\gamma}, \eqn{\beta},
#' \eqn{\psi}, \eqn{\tau}, and log-density chains.
#'
#' @param res An object of class \code{bayesmr_het_fit}.
#' @param include.burnin Logical; if \code{TRUE} the burn-in iterations
#'   (when stored) are included. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} prints warnings during conversion.
#'   Default is \code{TRUE}.
#'
#' @return A \code{\link[coda]{mcmc}} object.
#'
#' @seealso \code{\link{bayesmr_het}()}, \code{\link{bayesmr_het_fit-class}},
#'   \code{\link[coda]{mcmc}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' res <- bayesmr_het(dat, control = bayesmr_control(nsim = 500, burnin = 500))
#' mc <- bayesmr_het_fit_to_mcmc(res@results[[1]])
#' plot(mc)
#' }
#'
#' @export
bayesmr_het_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@dim[["n"]]

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

#' Conversion of a \code{bayesmr_het_fit_list} Object to a \code{list}
#'
#' Converts a multi-chain \code{\link{bayesmr_het_fit_list-class}} object to a
#' plain list of matrices (one per chain).
#'
#' @param res An object of class \code{bayesmr_het_fit_list}.
#' @param include.burnin Logical; if \code{TRUE} the burn-in iterations
#'   (when stored) are included. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} prints warnings during conversion.
#'   Default is \code{TRUE}.
#'
#' @return A list of numeric matrices, one per chain, with columns
#'   \code{gamma}, \code{beta}, \code{psi}, \code{tau}, \code{loglik},
#'   \code{logprior}, \code{logpost}.
#'
#' @seealso \code{\link{bayesmr_het}()},
#'   \code{\link{bayesmr_het_fit_list-class}},
#'   \code{\link{bayesmr_het_fit_list_to_mcmc.list}()}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' res <- bayesmr_het(dat, control = bayesmr_control(nsim = 500, burnin = 500,
#'                                                   nchains = 2))
#' lst <- bayesmr_het_fit_list_to_list(res)
#' bayesplot::mcmc_trace(lst, pars = "gamma")
#' }
#'
#' @export
bayesmr_het_fit_list_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@results[[1]]@dim[["n"]]

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

#' Conversion of a \code{bayesmr_het_fit_list} Object to \code{mcmc.list}
#'
#' Converts a multi-chain \code{\link{bayesmr_het_fit_list-class}} object to a
#' \code{\link[coda]{mcmc.list}} for convergence diagnostics.
#'
#' @param res An object of class \code{bayesmr_het_fit_list}.
#' @param include.burnin Logical; if \code{TRUE} the burn-in iterations
#'   (when stored) are included. Default is \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} prints warnings during conversion.
#'   Default is \code{TRUE}.
#'
#' @return A \code{\link[coda]{mcmc.list}} object.
#'
#' @seealso \code{\link{bayesmr_het}()},
#'   \code{\link{bayesmr_het_fit_list-class}},
#'   \code{\link[coda]{mcmc.list}},
#'   \code{\link{bayesmr_het_fit_list_to_list}()}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' res <- bayesmr_het(dat, control = bayesmr_control(nsim = 500, burnin = 500,
#'                                                   nchains = 2))
#' mc <- bayesmr_het_fit_list_to_mcmc.list(res)
#' coda::gelman.diag(mc)
#' }
#'
#' @export
bayesmr_het_fit_list_to_mcmc.list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  out <- bayesmr_het_fit_list_to_list(res, include.burnin = include.burnin, verbose = verbose)
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

#' Conversion of a \code{bayesmr_mix_fit} object to an object of class \code{mcmc}.
#'
#' \code{bayesmr_mix_fit_to_mcmc} converts an object of class \code{bayesmr_mix_fit}
#'   to a \code{coda::mcmc} object containing the scalar parameter chains
#'   (\eqn{\gamma}, \eqn{\alpha}) and log-density values. The per-SNP
#'   \eqn{\beta} and \eqn{\xi} chains are not included; use the dedicated
#'   post-processing functions in \code{bnpmr_postprocess.R} for those.
#'
#' @param res An object of type \code{bayesmr_mix_fit}.
#' @param include.burnin A logical scalar.
#' @param verbose A logical scalar.
#' @return An object of type \code{mcmc}.
#' @seealso \code{\link{bayesmr_mix}()}, \code{\link{bayesmr_mix_fit-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @export
bayesmr_mix_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

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
  theta[, 2] <- res@alpha.chain[tokeep]
  theta_nm[2] <- "alpha"
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

#' Conversion of a \code{bayesmr_mix_fit_list} object to a \code{list}.
#'
#' @param res An object of type \code{bayesmr_mix_fit_list}.
#' @param include.burnin A logical scalar.
#' @param verbose A logical scalar.
#' @return A list of matrices (one per chain).
#' @seealso \code{\link{bayesmr_mix}()}, \code{\link{bayesmr_mix_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @export
bayesmr_mix_fit_list_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

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
    theta <- matrix(NA, nrow = length(tokeep), ncol = 2 + 3)
    theta_nm <- character(2 + 3)
    theta[, 1] <- res@results[[c]]@gamma.chain[tokeep]
    theta_nm[1] <- "gamma"
    theta[, 2] <- res@results[[c]]@alpha.chain[tokeep]
    theta_nm[2] <- "alpha"
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

#' Conversion of a \code{bayesmr_mix_fit_list} object to an object of class \code{mcmc.list}.
#'
#' @param res An object of type \code{bayesmr_mix_fit_list}.
#' @param include.burnin A logical scalar.
#' @param verbose A logical scalar.
#' @return An object of type \code{mcmc.list}.
#' @seealso \code{\link{bayesmr_mix}()}, \code{\link{bayesmr_mix_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @export
bayesmr_mix_fit_list_to_mcmc.list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  out <- bayesmr_mix_fit_list_to_list(res, include.burnin = include.burnin, verbose = verbose)
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

#' Conversion of a \code{bayesmr_mix_het_fit} object to an object of class \code{mcmc}.
#'
#' \code{bayesmr_mix_het_fit_to_mcmc} converts an object of class
#'   \code{bayesmr_mix_het_fit} to a \code{coda::mcmc} object containing
#'   the scalar parameter chains (\eqn{\gamma}, \eqn{\alpha}, \eqn{\psi},
#'   \eqn{\tau}) and log-density values.
#'
#' @param res An object of type \code{bayesmr_mix_het_fit}.
#' @param include.burnin A logical scalar.
#' @param verbose A logical scalar.
#' @return An object of type \code{mcmc}.
#' @seealso \code{\link{bayesmr_mix_het}()}, \code{\link{bayesmr_mix_het_fit-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @export
bayesmr_mix_het_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

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
  theta[, 2] <- res@alpha.chain[tokeep]
  theta_nm[2] <- "alpha"
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

#' Conversion of a \code{bayesmr_mix_het_fit_list} object to a \code{list}.
#'
#' @param res An object of type \code{bayesmr_mix_het_fit_list}.
#' @param include.burnin A logical scalar.
#' @param verbose A logical scalar.
#' @return A list of matrices (one per chain).
#' @seealso \code{\link{bayesmr_mix_het}()}, \code{\link{bayesmr_mix_het_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @export
bayesmr_mix_het_fit_list_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

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
    theta <- matrix(NA, nrow = length(tokeep), ncol = 4 + 3)
    theta_nm <- character(4 + 3)
    theta[, 1] <- res@results[[c]]@gamma.chain[tokeep]
    theta_nm[1] <- "gamma"
    theta[, 2] <- res@results[[c]]@alpha.chain[tokeep]
    theta_nm[2] <- "alpha"
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

#' Conversion of a \code{bayesmr_mix_het_fit_list} object to an object of class \code{mcmc.list}.
#'
#' @param res An object of type \code{bayesmr_mix_het_fit_list}.
#' @param include.burnin A logical scalar.
#' @param verbose A logical scalar.
#' @return An object of type \code{mcmc.list}.
#' @seealso \code{\link{bayesmr_mix_het}()}, \code{\link{bayesmr_mix_het_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @export
bayesmr_mix_het_fit_list_to_mcmc.list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  out <- bayesmr_mix_het_fit_list_to_list(res, include.burnin = include.burnin, verbose = verbose)
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
