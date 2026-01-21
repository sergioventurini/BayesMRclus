#' Fitting function for BayesMR models with fixed heterogeneity.
#'
#' \code{bayesmr_fit()} is the main function that estimates a BayesMR model.
#'
#' @param data An object of class \code{\link{bayesmr_data}}).
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
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
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
bayesmr_fit <- function(data, p, G, control, prior, start) {
	n <- data@n
  data_obs <- data@data
	totiter <- control[["burnin"]] + control[["nsim"]]
	p <- as.integer(p)
	G <- as.integer(G)
	
	gamma.chain <- beta.chain <- array(NA, dim = totiter)
	loglik <- logprior <- logpost <- numeric(totiter)
	
	# recover prior hyperparameters
  hyper.gammaj.psi2 <- prior[["gammaj"]][["psi2"]]
  hyper.Gammaj.tau2 <- prior[["Gammaj"]][["tau2"]]
  hyper.gamma.mean <- prior[["gamma"]][["mean"]]
  hyper.gamma.var <- prior[["gamma"]][["var"]]
  hyper.beta.mean <- prior[["beta"]][["mean"]]
  hyper.beta.var <- prior[["beta"]][["var"]]
	
	# start iteration
	if (control[["verbose"]]) message("Running the MCMC simulation...")
	
	res.mcmc <- .Call('bayesmr_mcmc', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
		radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
		rn = as.integer(n),
		rp = as.integer(p),
		rG = as.integer(G),
		rtotiter = as.integer(totiter),
		rC_beta = as.double(control[["beta.prop"]]),
    rhyper_gammaj_psi2 = as.double(hyper.gammaj.psi2),
    rhyper_Gammaj_tau2 = as.double(hyper.Gammaj.tau2),
    rhyper_gamma_mean = as.double(hyper.gamma.mean),
    rhyper_gamma_var = as.double(hyper.gamma.var),
    rhyper_beta_mean = as.double(hyper.beta.mean),
    rhyper_beta_var = as.double(hyper.beta.var),
		rverbose = as.integer(control[["verbose"]])
	)

	gamma.chain <- array(res.mcmc[[1]], totiter)
  beta.chain <- array(res.mcmc[[2]], totiter)
	accept <- as.numeric(res.mcmc[[3]])
	loglik <- as.numeric(res.mcmc[[4]])
	logprior <- as.numeric(res.mcmc[[5]])
	logpost <- as.numeric(res.mcmc[[6]])

  # apply thinning
  if (control[["thin"]] > 1) {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = control[["thin"]])
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = control[["thin"]])
    }
  } else {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = 1)
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = 1)
    }
  }
  gamma.chain <- gamma.chain[tokeep, drop = FALSE]
  beta.chain <- beta.chain[tokeep, drop = FALSE]
  loglik <- loglik[tokeep]
  logprior <- logprior[tokeep]
  logpost <- logpost[tokeep]

  # return results
	out <- new("bayesmr_fit",
		gamma.chain = gamma.chain,
    beta.chain = beta.chain,
		accept = accept,
		data = data_obs,
		dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
    control = control,
    prior = prior,
    dim = list(n = n, p = p, G = G),
    model = new("bayesmr_model", p = p, G = G)
	)

	return(out)
}

#' Fitting function for BayesMR models with random heterogeneity.
#'
#' \code{bayesmr_het_fit()} is the main function that estimates a BayesMR model.
#'
#' @param data An object of class \code{\link{bayesmr_data}}).
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
#' @return A \code{bayesmr_het_fit_list} object.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bayesmr_data}} for a description of the data format.
#' @seealso \code{\link{bayesmr_fit_list}} for a description of the elements
#'   included in the returned object.
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
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
#' sim.bayesmr <- bayesmr_het(simdiss, p, G, control)
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
bayesmr_het_fit <- function(data, p, G, control, prior, start) {
  n <- data@n
  data_obs <- data@data
  totiter <- control[["burnin"]] + control[["nsim"]]
  p <- as.integer(p)
  G <- as.integer(G)
  
  gamma.chain <- beta.chain <- psi.chain <- tau.chain <- array(NA, dim = totiter)
  loglik <- logprior <- logpost <- numeric(totiter)
  
  # recover prior hyperparameters
  hyper.psi.alpha <- prior[["psi"]][["alpha"]]
  hyper.psi.nu <- prior[["psi"]][["nu"]]
  hyper.tau.alpha <- prior[["tau"]][["alpha"]]
  hyper.tau.nu <- prior[["tau"]][["nu"]]
  hyper.gamma.mean <- prior[["gamma"]][["mean"]]
  hyper.gamma.var <- prior[["gamma"]][["var"]]
  hyper.beta.mean <- prior[["beta"]][["mean"]]
  hyper.beta.var <- prior[["beta"]][["var"]]
  
  # start iteration
  if (control[["verbose"]]) message("Running the MCMC simulation...")
  
  res.mcmc <- .Call('bayesmr_mcmc_het', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
    radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
    radpsi = as.double(start$psi),
    radtau = as.double(start$tau),
    rn = as.integer(n),
    rp = as.integer(p),
    rG = as.integer(G),
    rtotiter = as.integer(totiter),
    rC_beta = as.double(control[["beta.prop"]]),
    rC_psi = as.double(control[["psi.prop"]]),
    rC_tau = as.double(control[["tau.prop"]]),
    rhyper_gamma_mean = as.double(hyper.gamma.mean),
    rhyper_gamma_var = as.double(hyper.gamma.var),
    rhyper_beta_mean = as.double(hyper.beta.mean),
    rhyper_beta_var = as.double(hyper.beta.var),
    rhyper_psi_alpha = as.double(hyper.psi.alpha),
    rhyper_psi_nu = as.double(hyper.psi.nu),
    rhyper_tau_alpha = as.double(hyper.tau.alpha),
    rhyper_tau_nu = as.double(hyper.tau.nu),
    rverbose = as.integer(control[["verbose"]])
  )

  gamma.chain <- array(res.mcmc[[1]], totiter)
  beta.chain <- array(res.mcmc[[2]], totiter)
  psi.chain <- array(res.mcmc[[3]], totiter)
  tau.chain <- array(res.mcmc[[4]], totiter)
  accept <- t(array(res.mcmc[[5]], c(G, 2)))  # CHECK THIS!!!
  loglik <- as.numeric(res.mcmc[[6]])
  logprior <- as.numeric(res.mcmc[[7]])
  logpost <- as.numeric(res.mcmc[[8]])

  # apply thinning
  if (control[["thin"]] > 1) {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = control[["thin"]])
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = control[["thin"]])
    }
  } else {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = 1)
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = 1)
    }
  }
  gamma.chain <- gamma.chain[tokeep, drop = FALSE]
  beta.chain <- beta.chain[tokeep, drop = FALSE]
  psi.chain <- psi.chain[tokeep, drop = FALSE]
  tau.chain <- tau.chain[tokeep, drop = FALSE]
  loglik <- loglik[tokeep]
  logprior <- logprior[tokeep]
  logpost <- logpost[tokeep]

  # return results
  out <- new("bayesmr_het_fit",
    gamma.chain = gamma.chain,
    beta.chain = beta.chain,
    psi.chain = psi.chain,
    tau.chain = tau.chain,
    accept = accept,
    data = data_obs,
    dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
    control = control,
    prior = prior,
    dim = list(n = n, p = p, G = G),
    model = new("bayesmr_model", p = p, G = G)
  )

  return(out)
}

#' Fitting function for BayesMR models with fixed heterogeneity.
#'
#' \code{bayesmr_fit()} is the main function that estimates a BayesMR model.
#'
#' @param data An object of class \code{\link{bayesmr_data}}).
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
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
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
bayesmr_mix_fit <- function(data, p, G, control, prior, start) {
  n <- data@n
  data_obs <- data@data
  totiter <- control[["burnin"]] + control[["nsim"]]
  p <- as.integer(p)
  G <- as.integer(G)
  
  gamma.chain <- alpha.chain <- array(NA, dim = totiter)
  beta.chain <- xi.chain <- array(NA, dim = c(totiter, n))
  loglik <- logprior <- logpost <- numeric(totiter)
  
  # recover prior hyperparameters
  hyper.gammaj.psi2 <- prior[["gammaj"]][["psi2"]]
  hyper.Gammaj.tau2 <- prior[["Gammaj"]][["tau2"]]
  hyper.gamma.mean <- prior[["gamma"]][["mean"]]
  hyper.gamma.var <- prior[["gamma"]][["var"]]
  hyper.beta.mean <- prior[["beta"]][["mean"]]
  hyper.beta.var <- prior[["beta"]][["var"]]
  hyper.alpha.a <- prior[["alpha"]][["a"]]
  hyper.alpha.b <- prior[["alpha"]][["b"]]
  
  # start iteration
  if (control[["verbose"]]) message("Running the MCMC simulation...")
  
  res.mcmc <- .Call('bayesmr_mix_mcmc', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
    radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
    radxi = as.integer(start$xi),
    radalpha = as.double(start$alpha),
    radK = as.double(start$K),
    rn = as.integer(n),
    rp = as.integer(p),
    rG = as.integer(G),
    rtotiter = as.integer(totiter),
    rC_beta = as.double(control[["beta.prop"]]),
    rm_beta = as.double(control[["beta.m"]]),
    rhyper_gammaj_psi2 = as.double(hyper.gammaj.psi2),
    rhyper_Gammaj_tau2 = as.double(hyper.Gammaj.tau2),
    rhyper_gamma_mean = as.double(hyper.gamma.mean),
    rhyper_gamma_var = as.double(hyper.gamma.var),
    rhyper_beta_mean = as.double(hyper.beta.mean),
    rhyper_beta_var = as.double(hyper.beta.var),
    rhyper_alpha_a = as.double(hyper.alpha.a),
    rhyper_alpha_b = as.double(hyper.alpha.b),
    rverbose = as.integer(control[["verbose"]])
  )

  gamma.chain <- array(res.mcmc[[1]], totiter)
  beta.chain <- t(array(res.mcmc[[2]], c(n, totiter)))
  xi.chain <- t(array(res.mcmc[[3]], c(n, totiter)))
  alpha.chain <- array(res.mcmc[[4]], totiter)
  accept <- as.numeric(res.mcmc[[5]])
  loglik <- as.numeric(res.mcmc[[6]])
  logprior <- as.numeric(res.mcmc[[7]])
  logpost <- as.numeric(res.mcmc[[8]])

  # apply thinning
  if (control[["thin"]] > 1) {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = control[["thin"]])
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = control[["thin"]])
    }
  } else {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = 1)
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = 1)
    }
  }
  gamma.chain <- gamma.chain[tokeep, drop = FALSE]
  beta.chain <- beta.chain[tokeep, , drop = FALSE]
  xi.chain <- xi.chain[tokeep, , drop = FALSE]
  alpha.chain <- alpha.chain[tokeep, drop = FALSE]
  loglik <- loglik[tokeep]
  logprior <- logprior[tokeep]
  logpost <- logpost[tokeep]

  # return results
  out <- new("bayesmr_mix_fit",
    gamma.chain = gamma.chain,
    beta.chain = beta.chain,
    xi.chain = xi.chain,
    alpha.chain = alpha.chain,
    accept = accept,
    data = data_obs,
    dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
    control = control,
    prior = prior,
    dim = list(n = n, p = p, G = G),
    model = new("bayesmr_model", p = p, G = G)
  )

  return(out)
}
