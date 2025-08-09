library(BayesMRclus)

# prepare data
data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])
n <- nrow(data)
zhaodata <- new("bayesmr_data", data = data, n = n)
summary(zhaodata)
plot(zhaodata)

beta_z <- 0.3
beta_se_z <- 0.2
tau2_z <- 1e-3
hpar <- list(mu_gamma = mean(data[, 1]),
             sigma2_gamma = var(data[, 1])/nrow(data),
             mu_beta = beta_z,
             sigma2_beta = beta_se_z^2,
             psi2 = mean(data[, 3]^2),
             tau2 = tau2_z)

iter <- burnin + nsim
start <- rnorm(2)
tune <- .1

seed <- 2301
set.seed(seed)

res_R <- mcmc_bayesmr(data, hpar, iter, start, tune)
summary(res_R$draws[, 1])
plot(res_R$draws[, 1], type = "l")
summary(res_R$draws[, 2])
plot(res_R$draws[, 2], type = "l")

###

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 1e-2),
                       Gammaj = list(tau2 = 1e-3),
                       gamma = list(mean = 0,
                                    var = 1e-1),
                       beta = list(mean = 0, var = 1e1))

gamma <- seq(-.05, .05, .01)
beta <- seq(1, 4, .01)
out <- gamma_beta_post(gamma, beta, data, prior, log = FALSE, verbose = FALSE)
str(out)
contour(gamma, beta, out, xlab = "gamma", ylab = "beta")
abline(v = prior$gamma$mean, lty = 2)
abline(h = prior$beta$mean, lty = 2)
summary(as.vector(out))
# persp(gamma, beta, out)

###

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 1e0),
                       Gammaj = list(tau2 = 1e-3),
                       gamma = list(mean = 0,
                                    var = 1e-1),
                       beta = list(mean = 1, var = 1e1))

gamma <- 0
beta <- seq(-1e1, 1e1, .01)
out <- gamma_beta_post(gamma, beta, data, prior, log = FALSE, verbose = FALSE)
str(out)
plot(beta, out, xlab = "beta", ylab = "full conditional distribution",
     type = "l")
abline(v = prior$beta$mean, lty = 2)
summary(as.vector(out))

###

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 1e0),
                       Gammaj = list(tau2 = 1e-3),
                       gamma = list(mean = 0,
                                    var = 1e-1),
                       beta = list(mean = 1, var = 1e1))

gamma <- 0
beta_min <- -1.5
beta_max <- 1.5
beta_step <- 0.001
beta_seq <- seq(beta_min, beta_max, beta_step)
out <- beta_marg_post_drv(beta_seq, gamma, data, prior)
plot(beta_seq, out, type = "l")
abline(h = 0, lty = 1, col = "gray")
abline(v = prior$beta$mean, lty = 2, col = "gray")

# uniroot(beta_marg_post_drv, interval = range(beta_seq),
#         gamma = gamma, data = data, prior = prior,
#         tol = 1e-5, maxiter = 1000, trace = 2)
roots <- find_all_roots(beta_marg_post_drv,
                        lower = beta_min, upper = beta_max,
                        gamma = gamma, data = data, prior = prior,
                        n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                        eps_small = 1e-8)
roots
abline(v = roots, lty = 2, col = "gray")
