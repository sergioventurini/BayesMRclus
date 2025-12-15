# needed because optim() requires a function returning a scalar
bayesmr_logLik_optim <- function(par, data, prior) {
  gamma <- par[1]
  beta  <- par[2]

  as.numeric(
    bayesmr_logLik(
      gamma = gamma,
      beta  = beta,
      data  = data,
      prior = prior,
      log   = TRUE
    )
  )
}

# needed because optim() requires a function returning a scalar
gamma_beta_logpost_optim <- function(par, data, prior) {
  gamma <- par[1]
  beta  <- par[2]

  as.numeric(
    gamma_beta_post(
      gamma = gamma,
      beta  = beta,
      data  = data,
      prior = prior,
      log   = TRUE
    )
  )
}

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
bayesmr_noclus_mle <- function(
  data,
  prior,
  init = c(gamma = 0, beta = 0),
  method = "BFGS",
  control = list(fnscale = -1, reltol = 1e-10, maxit = 1000)) {
  opt <- optim(
    par = init,
    fn = bayesmr_logLik_optim,
    data = data,
    prior = prior,
    method = method,
    hessian = TRUE,
    control = control
  )

  vcov <- tryCatch(
    solve(-opt$hessian),
    error = function(e) matrix(NA, 2, 2)
  )

  # return results
  list(
    par = setNames(opt$par, c("gamma", "beta")),
    value = opt$value,
    convergence = opt$convergence,
    hessian = opt$hessian,
    vcov = vcov,
    sd = sqrt(diag(vcov)),
    optim = opt
  )
}

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
bayesmr_noclus_optim <- function(
  data,
  prior,
  init = c(gamma = 0, beta = 0),
  method = "BFGS",
  control = list(fnscale = -1, reltol = 1e-10, maxit = 1000)) {
  opt <- optim(
    par = init,
    fn = gamma_beta_logpost_optim,
    data = data,
    prior = prior,
    method = method,
    hessian = TRUE,
    control = control
  )

  vcov <- tryCatch(
    solve(-opt$hessian),
    error = function(e) matrix(NA, 2, 2)
  )

  # return results
  list(
    par = setNames(opt$par, c("gamma", "beta")),
    value = opt$value,
    convergence = opt$convergence,
    hessian = opt$hessian,
    vcov = vcov,
    sd = sqrt(diag(vcov)),
    optim = opt
  )
}
