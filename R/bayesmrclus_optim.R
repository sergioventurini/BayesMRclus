#' Alternating conditional maximization of the bivariate posterior for a BayesMR model with no clustering.
#'
#' \code{bayesmr_fit()} is the main function that estimates a BayesMR model.
#'
#' @param data An object of class \code{\link{bayesmr_data}})..
#' @param p A length-one numeric vector indicating the number of dimensions of the
#'   latent space.
#' @param G A length-one numeric vector indicating the number of cluster to
#'   partition the \emph{S} subjects.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution See
#'   \code{\link{bayesmr_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{bayesmr_prior}()} for more details.
#' @param start A named list of starting values for the MCMC algorithm (see
#'   \code{\link{bayesmr_init}}).
#' @return A \code{bayesmr_fit_list} object.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bayesmr_data}} for a description of the data format.
#' @seealso \code{\link{bayesmr_fit_list}} for a description of the elements
#'   included in the returned object.
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2024), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 20000
#' nsim <- 10000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
#'
#' summary(sim.bayesmr, include.burnin = FALSE)
#'
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("teal")
#' plot(sim.bayesmr, what = "trace", regex_pars = "eta")
#'
#' z <- bayesmr_get_configuration(sim.bayesmr, chain = 1, est = "mean",
#'   labels = 1:16)
#' summary(z)
#' color_scheme_set("mix-pink-blue")
#' graph <- plot(z, size = 2, size_lbl = 3)
#' graph + panel_bg(fill = "gray90", color = NA)
#' }
#' @export
bayesmr_noclus_optim <- function(data, prior, start = rep(0, 2), maxiter = 1000, tol = 1e-10,
                                 beta_min = -20, beta_max = 20, beta_step = 0.001,
                                 n = 1000, tol_x = 1e-10, tol_f = 1e-12, eps_small = 1e-8,
                                 verbose = TRUE) {
	gamma_chain <- beta_chain <- logpost <- NA

  # start iterations
  if (verbose) message("Running the optimization procedure...")
  
  gamma_chain <- start[1]  # gamma starting value
  beta_chain <- start[2]   # beta starting value
  logpost <- gamma_beta_post(gamma_chain[1], beta_chain[1], data, prior, log = TRUE, verbose = FALSE)
  if (verbose)
    message("  - iteration ", 1, " - gamma value: ", round(gamma_chain[1], digits = 6),
            " - beta value: ", round(beta_chain[1], digits = 6),
            " - log posterior value: ", format(logpost, nsamll = 8, scientific = FALSE))
  for (i in 2:maxiter) {
    # find gamma mode conditionally on beta
    gamma_optim <- gamma_marg_post_mode(beta_chain[i - 1], data, prior)
    gamma_chain <- c(gamma_chain, gamma_optim$modes)

    # find beta mode(s) conditionally on gamma
    beta_optim <- beta_marg_post_mode(gamma_chain[i], data, prior,
                                      beta_min = beta_min, beta_max = beta_max,
                                      beta_step = beta_step, n = n, tol_x = tol_x,
                                      tol_f = tol_f, eps_small = eps_small, log = TRUE)
    if (length(beta_optim$modes) == 1) {
      beta_chain <- c(beta_chain, beta_optim$modes)
    }
    else {
      if (identical(beta_optim$dens[1], beta_optim$dens[2]))
        stop("the joint posterior has two modes.")
      beta_chain <- c(beta_chain, beta_optim$modes[which.max(beta_optim$dens)])
    }

    logdens <- gamma_beta_post(gamma_chain[i], beta_chain[i], data, prior, log = TRUE, verbose = FALSE)
    logpost <- c(logpost, logdens)

    if (verbose) message("  - iteration ", i, " - gamma value: ", round(gamma_chain[i], digits = 6),
                         " - beta value: ", round(beta_chain[i], digits = 6),
                         " - log posterior value: ", format(logdens, nsamll = 8, scientific = FALSE))

    # check convergence in terms of change in parameter values
    if (max(abs(diff(gamma_chain[(i - 1):i])), abs(diff(beta_chain[(i - 1):i]))) < tol) {
      break
    }
  }

  # return results
	out <- list(gamma_best = gamma_chain[i], beta_best = beta_chain[i],
              gamma_chain = gamma_chain, beta_chain = beta_chain,
              logpost = logpost, iter = i)

	out
}
