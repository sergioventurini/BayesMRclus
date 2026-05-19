#' Fit a BayesMR model with fixed heterogeneity.
#'
#' `bayesmr_fit()` runs the MCMC sampler for a BayesMR model with fixed
#' heterogeneity and returns the fitted chain. This is a lower-level fitting
#' function called by [`bayesmr()`]; most users should call [`bayesmr()`]
#' directly.
#'
#' @param data An object of class [`bayesmr_data-class`] containing the summary
#'   data to be analyzed.
#' @param control A named list of control parameters for the MCMC sampler, such
#'   as the number of iterations, burn-in, thinning, proposal variance, and
#'   verbosity. See [`bayesmr_control()`] for details.
#' @param prior A named list containing the prior hyperparameters. See
#'   [`bayesmr_prior()`] for details.
#' @param start A named list of starting values for the MCMC algorithm, typically
#'   created by [`bayesmr_init()`]. It must contain at least the elements
#'   `gamma` and `beta`.
#'
#' @return An object of class [`bayesmr_fit-class`] containing the simulated
#'   chains, acceptance information, log-density values, data, control settings,
#'   prior specification, and model dimensions.
#'
#' @details
#' The function calls the compiled sampler `bayesmr_mcmc_noclus_wrap`, then
#' applies burn-in removal and thinning according to the supplied `control`
#' settings. If `control$store.burnin` is `TRUE`, burn-in iterations are retained
#' before thinning; otherwise only post-burn-in draws are stored.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso
#' [`bayesmr()`] for the main user-facing fitting function;
#' [`bayesmr_data-class`] for the input data format;
#' [`bayesmr_fit-class`] for the returned object;
#' [`bayesmr_control()`], [`bayesmr_prior()`], and [`bayesmr_init()`].
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#'
#' ctrl <- bayesmr_control(nsim = 500, burnin = 500)
#' pr <- bayesmr_prior()
#' st <- bayesmr_init(
#'   data = dat,
#'   random_start = ctrl$random_start,
#'   K_start = ctrl$K_start
#' )
#'
#' fit <- bayesmr_fit(
#'   data = dat,
#'   control = ctrl,
#'   prior = pr,
#'   start = st
#' )
#' }
#'
#' @export
bayesmr_fit <- function(data, control, prior, start) {
	n <- data@n
  data_obs <- data@data
	totiter <- control[["burnin"]] + control[["nsim"]]
	
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
	
	res.mcmc <- .Call('bayesmr_mcmc_noclus_wrap', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
		radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
		rn = as.integer(n),
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
    dim = list(n = n)
	)

	return(out)
}

#' Fit a BayesMR model with random heterogeneity.
#'
#' `bayesmr_het_fit()` runs the MCMC sampler for a BayesMR model with random
#' heterogeneity and returns the fitted chain. This is a lower-level fitting
#' function; most users should call [`bayesmr()`] or the relevant user-facing
#' wrapper directly.
#'
#' @param data An object of class [`bayesmr_data-class`] containing the summary
#'   data to be analyzed.
#' @param control A named list of control parameters for the MCMC sampler, such
#'   as the number of iterations, burn-in, thinning, proposal variances, and
#'   verbosity. See [`bayesmr_control()`] for details.
#' @param prior A named list containing the prior hyperparameters. See
#'   [`bayesmr_prior()`] for details.
#' @param start A named list of starting values for the MCMC algorithm, typically
#'   created by [`bayesmr_init()`]. It must contain at least the elements
#'   `gamma`, `beta`, `psi`, and `tau`.
#'
#' @return An object of class [`bayesmr_het_fit-class`] containing the simulated
#'   chains for `gamma`, `beta`, `psi`, and `tau`, acceptance information,
#'   log-density values, data, control settings, prior specification, and model
#'   dimensions.
#'
#' @details
#' The function calls the compiled sampler `bayesmr_mcmc_noclus_het_wrap`, then
#' applies burn-in removal and thinning according to the supplied `control`
#' settings. If `control$store.burnin` is `TRUE`, burn-in iterations are retained
#' before thinning; otherwise only post-burn-in draws are stored.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso
#' [`bayesmr()`] for the main user-facing fitting function;
#' [`bayesmr_data-class`] for the input data format;
#' [`bayesmr_het_fit-class`] for the returned object;
#' [`bayesmr_control()`], [`bayesmr_prior()`], and [`bayesmr_init()`].
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#'
#' ctrl <- bayesmr_control(nsim = 500, burnin = 500)
#' pr <- bayesmr_prior()
#' st <- bayesmr_init(
#'   data = dat,
#'   random_start = ctrl$random_start,
#'   K_start = ctrl$K_start
#' )
#'
#' fit <- bayesmr_het_fit(
#'   data = dat,
#'   control = ctrl,
#'   prior = pr,
#'   start = st
#' )
#' }
#'
#' @export
bayesmr_het_fit <- function(data, control, prior, start) {
  n <- data@n
  data_obs <- data@data
  totiter <- control[["burnin"]] + control[["nsim"]]
  
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
  
  res.mcmc <- .Call('bayesmr_mcmc_noclus_het_wrap', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
    radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
    radpsi = as.double(start$psi),
    radtau = as.double(start$tau),
    rn = as.integer(n),
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
  accept <- t(array(res.mcmc[[5]], 2))
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
    dim = list(n = n)
  )

  return(out)
}

#' Fit a BayesMR mixture model.
#'
#' `bayesmr_mix_fit()` runs the MCMC sampler for a BayesMR mixture model with
#' cluster-specific causal effects and fixed heterogeneity. This is a lower-level
#' fitting function; most users should call [`bayesmr()`] or the relevant
#' user-facing wrapper directly.
#'
#' @param data An object of class [`bayesmr_data-class`] containing the summary
#'   data to be analyzed.
#' @param control A named list of control parameters for the MCMC sampler, such
#'   as the number of iterations, burn-in, thinning, proposal variance, mixture
#'   proposal settings, and verbosity. See [`bayesmr_control()`] for details.
#' @param prior A named list containing the prior hyperparameters. See
#'   [`bayesmr_prior()`] for details.
#' @param start A named list of starting values for the MCMC algorithm, typically
#'   created by [`bayesmr_init()`]. It must contain at least the elements
#'   `gamma`, `beta`, `xi`, `alpha`, and `K`.
#'
#' @return An object of class [`bayesmr_mix_fit-class`] containing the simulated
#'   chains for `gamma`, `beta`, `xi`, and `alpha`, acceptance information,
#'   log-density values, data, control settings, prior specification, and model
#'   dimensions.
#'
#' @details
#' The function calls the compiled sampler `bayesmr_mcmc_mix_wrap`, then applies
#' burn-in removal and thinning according to the supplied `control` settings. If
#' `control$store.burnin` is `TRUE`, burn-in iterations are retained before
#' thinning; otherwise only post-burn-in draws are stored.
#'
#' The returned `beta.chain` and `xi.chain` matrices have one row per stored MCMC
#' iteration and one column per genetic instrument.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso
#' [`bayesmr()`] for the main user-facing fitting function;
#' [`bayesmr_data-class`] for the input data format;
#' [`bayesmr_mix_fit-class`] for the returned object;
#' [`bayesmr_control()`], [`bayesmr_prior()`], and [`bayesmr_init()`].
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#'
#' ctrl <- bayesmr_control(nsim = 500, burnin = 500)
#' pr <- bayesmr_prior()
#' st <- bayesmr_init(
#'   data = dat,
#'   random_start = ctrl$random_start,
#'   K_start = ctrl$K_start
#' )
#'
#' fit <- bayesmr_mix_fit(
#'   data = dat,
#'   control = ctrl,
#'   prior = pr,
#'   start = st
#' )
#' }
#'
#' @export
bayesmr_mix_fit <- function(data, control, prior, start) {
  n <- data@n
  data_obs <- data@data
  totiter <- control[["burnin"]] + control[["nsim"]]
  
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
  
  res.mcmc <- .Call('bayesmr_mcmc_mix_wrap', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
    radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
    radxi = as.integer(start$xi),
    radalpha = as.double(start$alpha),
    radK = as.double(start$K),
    rn = as.integer(n),
    rtotiter = as.integer(totiter),
    rC_beta = as.double(control[["beta.prop"]]),
    rm_beta = as.integer(control[["beta.m"]]),
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
    dim = list(n = n)
  )

  return(out)
}

#' Fit a BayesMR mixture model with random heterogeneity.
#'
#' `bayesmr_mix_het_fit()` runs the MCMC sampler for a BayesMR mixture model
#' with cluster-specific causal effects and random heterogeneity. This is a
#' lower-level fitting function; most users should call [`bayesmr()`] or the
#' relevant user-facing wrapper directly.
#'
#' @param data An object of class [`bayesmr_data-class`] containing the summary
#'   data to be analyzed.
#' @param control A named list of control parameters for the MCMC sampler, such
#'   as the number of iterations, burn-in, thinning, proposal variances, mixture
#'   proposal settings, and verbosity. See [`bayesmr_control()`] for details.
#' @param prior A named list containing the prior hyperparameters. See
#'   [`bayesmr_prior()`] for details.
#' @param start A named list of starting values for the MCMC algorithm, typically
#'   created by [`bayesmr_init()`]. It must contain at least the elements
#'   `gamma`, `beta`, `xi`, `alpha`, `psi`, `tau`, and `K`.
#'
#' @return An object of class [`bayesmr_mix_het_fit-class`] containing the
#'   simulated chains for `gamma`, `beta`, `xi`, `alpha`, `psi`, and `tau`,
#'   acceptance information, log-density values, data, control settings, prior
#'   specification, and model dimensions.
#'
#' @details
#' The function calls the compiled sampler `bayesmr_mcmc_mix_het_wrap`, then
#' applies burn-in removal and thinning according to the supplied `control`
#' settings. If `control$store.burnin` is `TRUE`, burn-in iterations are retained
#' before thinning; otherwise only post-burn-in draws are stored.
#'
#' The returned `beta.chain` and `xi.chain` matrices have one row per stored MCMC
#' iteration and one column per genetic instrument.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso
#' [`bayesmr()`] for the main user-facing fitting function;
#' [`bayesmr_data-class`] for the input data format;
#' [`bayesmr_mix_het_fit-class`] for the returned object;
#' [`bayesmr_control()`], [`bayesmr_prior()`], and [`bayesmr_init()`].
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#'
#' ctrl <- bayesmr_control(nsim = 500, burnin = 500)
#' pr <- bayesmr_prior()
#' st <- bayesmr_init(
#'   data = dat,
#'   random_start = ctrl$random_start,
#'   K_start = ctrl$K_start
#' )
#'
#' fit <- bayesmr_mix_het_fit(
#'   data = dat,
#'   control = ctrl,
#'   prior = pr,
#'   start = st
#' )
#' }
#'
#' @export
bayesmr_mix_het_fit <- function(data, control, prior, start) {
  n <- data@n
  data_obs <- data@data
  totiter <- control[["burnin"]] + control[["nsim"]]
  
  gamma.chain <- alpha.chain <- psi.chain <- tau.chain <- array(NA, dim = totiter)
  beta.chain <- xi.chain <- array(NA, dim = c(totiter, n))
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
  hyper.alpha.a <- prior[["alpha"]][["a"]]
  hyper.alpha.b <- prior[["alpha"]][["b"]]
  
  # start iteration
  if (control[["verbose"]]) message("Running the MCMC simulation...")
  
  res.mcmc <- .Call('bayesmr_mcmc_mix_het_wrap', PACKAGE = 'BayesMRclus',
    radData = as.double(unlist(data_obs)),
    radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
    radxi = as.integer(start$xi),
    radalpha = as.double(start$alpha),
    radpsi = as.double(start$psi),
    radtau = as.double(start$tau),
    radK = as.double(start$K),
    rn = as.integer(n),
    rtotiter = as.integer(totiter),
    rC_beta = as.double(control[["beta.prop"]]),
    rC_logpsi = as.double(control[["psi.prop"]]),
    rC_logtau = as.double(control[["tau.prop"]]),
    rm_beta = as.integer(control[["beta.m"]]),
    rhyper_gamma_mean = as.double(hyper.gamma.mean),
    rhyper_gamma_var = as.double(hyper.gamma.var),
    rhyper_beta_mean = as.double(hyper.beta.mean),
    rhyper_beta_var = as.double(hyper.beta.var),
    rhyper_psi_alpha = as.double(hyper.psi.alpha),
    rhyper_psi_nu = as.double(hyper.psi.nu),
    rhyper_tau_alpha = as.double(hyper.tau.alpha),
    rhyper_tau_nu = as.double(hyper.tau.nu),
    rhyper_alpha_a = as.double(hyper.alpha.a),
    rhyper_alpha_b = as.double(hyper.alpha.b),
    rverbose = as.integer(control[["verbose"]])
  )

  gamma.chain <- array(res.mcmc[[1]], totiter)
  beta.chain <- t(array(res.mcmc[[2]], c(n, totiter)))
  xi.chain <- t(array(res.mcmc[[3]], c(n, totiter)))
  alpha.chain <- array(res.mcmc[[4]], totiter)
  psi.chain <- array(res.mcmc[[5]], totiter)
  tau.chain <- array(res.mcmc[[6]], totiter)
  accept <- as.numeric(res.mcmc[[7]])
  loglik <- as.numeric(res.mcmc[[8]])
  logprior <- as.numeric(res.mcmc[[9]])
  logpost <- as.numeric(res.mcmc[[10]])

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
  psi.chain <- psi.chain[tokeep, drop = FALSE]
  tau.chain <- tau.chain[tokeep, drop = FALSE]
  loglik <- loglik[tokeep]
  logprior <- logprior[tokeep]
  logpost <- logpost[tokeep]

  # return results
  out <- new("bayesmr_mix_het_fit",
    gamma.chain = gamma.chain,
    beta.chain = beta.chain,
    xi.chain = xi.chain,
    alpha.chain = alpha.chain,
    psi.chain = psi.chain,
    tau.chain = tau.chain,
    accept = accept,
    data = data_obs,
    dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
    control = control,
    prior = prior,
    dim = list(n = n)
  )

  return(out)
}
