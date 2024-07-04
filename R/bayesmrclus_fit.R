#' Fitter function for BayesMR models.
#'
#' \code{bayesmr_fit()} is the main function that estimates a BayesMR model.
#'
#' @param D A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
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
#'   Venturini, S., Piccarreta, R. (2021), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{bayesmr}
#'   Package in \code{R}", Journal of Statistical Software, 100, 16, 1--35, <10.18637/jss.v100.i16>.
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
	n <- nrow(data)
	totiter <- control[["burnin"]] + control[["nsim"]]
	p <- as.integer(p)
	G <- as.integer(G)
	
	gamma.chain <- beta.chain <- array(NA, dim = totiter)
	loglik <- logprior <- logpost <- numeric(totiter)
	
	# recover prior hyperparameters
  hyper.gammaj.gamma <- prior[["gammaj"]][["gamma"]]
  hyper.gammaj.psi2 <- prior[["gammaj"]][["psi2"]]
  hyper.Gammaj.tau2 <- prior[["Gammaj"]][["tau2"]]
  hyper.gamma.mean <- prior[["gamma"]][["mean"]]
  hyper.gamma.var <- prior[["gamma"]][["var"]]
  hyper.beta.mean <- prior[["beta"]][["mean"]]
  hyper.beta.var <- prior[["beta"]][["var"]]
	
	# start iteration
	if (control[["verbose"]]) message("Running the MCMC simulation...")
	
	res.mcmc <- .Call('bayesmr_mcmc', PACKAGE = 'bayesmr',
    radData = as.double(unlist(data)),
		radgamma = as.double(start$gamma),
    radbeta = as.double(start$beta),
		rn = as.integer(n),
		rp = as.integer(p),
		rG = as.integer(G),
		rtotiter = as.integer(totiter),
		rsigma2_beta = as.double(control[["beta.prop"]]),
    rhyper_gammaj_gamma = as.double(hyper.gammaj.gamma),
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
	accept <- t(array(res.mcmc[[3]], c(G, 2)))
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
  gamma.chain <- gamma.chain[tokeep, , , , drop = FALSE]
  beta.chain <- beta.chain[tokeep, , , , drop = FALSE]
  loglik <- loglik[tokeep]
  logprior <- logprior[tokeep]
  logpost <- logpost[tokeep]

  # return results
	out <- new("bayesmr_fit",
		gamma.chain = gamma.chain,
    beta.chain = beta.chain,
		accept = accept,
		data = data,
		dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
    control = control,
    prior = prior,
    dim = list(n = n, p = p, G = G),
    model = new("bayesmr_model", p = p, G = G)
	)

	return(out)
}

#' Log-likelihood for BayesMR models.
#'
#' \code{bayesmr_logLik()} computes the log-likelihood value for a BayesMR model.
#'
#' @param data A data frame ...
#' @param gammaj A numeric vector ...
#' @param Gammaj A numeric vector ...
#'
#' @return A length-one numeric vector of the log-likelihood value.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso \code{\link{bayesmr}()}.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2021), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{bayesmr}
#'   Package in \code{R}", Journal of Statistical Software, 100, 16, 1--35, <10.18637/jss.v100.i16>.
#'
#' @export
bayesmr_logLik <- function(data, gammaj, Gammaj) {
	gammaj_hat <- data[, 1]
  Gammaj_hat <- data[, 2]
  sigmaj_X <- data[, 3]
  sigmaj_Y <- data[, 4]
  ll <- sum(dnorm(gammaj_hat, gammaj, sigmaj_X, log = TRUE)) +
        sum(dnorm(Gammaj_hat, Gammaj, sigmaj_Y, log = TRUE))

	return(ll)
}
