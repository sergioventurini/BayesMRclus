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
                       gamma = list(mean = 0, var = 1e-1),
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
                       gamma = list(mean = 0, var = 1e-1),
                       beta = list(mean = 1, var = 1e1))

gamma <- 0
beta <- seq(-1e1, 1e1, .01)
out <- gamma_beta_post(gamma, beta, data, prior, log = FALSE, verbose = FALSE)
str(out)
plot(beta, out, xlab = "beta", ylab = "full conditional distribution",
     type = "l")
abline(v = prior$beta$mean, lty = 2)
summary(as.vector(out))

### Find beta marginal posterior modes (given gamma)

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 1),
                       Gammaj = list(tau2 = 1),
                       gamma = list(mean = 0, var = 1),
                       beta = list(mean = 0, var = 1))

gamma <- 0
beta_min <- -5
beta_max <- 5
beta_step <- 0.001
beta_modes <- beta_marg_post_mode(gamma, data, prior,
                                  beta_min = beta_min, beta_max = beta_max,
                                  beta_step = beta_step,
                                  n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                                  eps_small = 1e-8)
beta_modes
beta_seq <- seq(beta_min, beta_max, beta_step)
out <- beta_marg_post_drv(beta_seq, gamma, data, prior)
plot(beta_seq, out, type = "l")
abline(h = 0, lty = 1, col = "gray")
abline(v = prior$beta$mean, lty = 2, col = "gray")
abline(v = beta_modes$modes, lty = 2, col = "gray")

### Find tau2 value where beta marginal posterior bifurcates (given gamma)

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 1e-1),
                       Gammaj = list(tau2 = 1e-2),
                       gamma = list(mean = 0, var = .1),
                       beta = list(mean = 0, var = .1))

tau2_min <- 0
tau2_max <- 5
tau2_step <- 0.0001
tau2_seq <- seq(tau2_min, tau2_max, tau2_step)

gamma <- 0
beta_min <- -5
beta_max <- 5
beta_step <- 0.001

res <- beta_mode_dist(gamma, data, prior,
                      tau2_min = tau2_min, tau2_max = tau2_max,
                      tau2_step = tau2_step,
                      beta_min = beta_min, beta_max = beta_max,
                      beta_step = beta_step,
                      n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                      eps_small = 1e-8, verbose = TRUE)

tau2_star <- res$tau2_star
iter <- res$iter
plot(tau2_seq[1:iter], res$dist,
     type = "l", xlab = "tau2", ylab = "distance",
     main = "Distance between beta marginal posterior modes")
iter
tau2_star

### Numerical optimization of the joint posterior

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 0.1),
                       Gammaj = list(tau2 = 0.1),
                       gamma = list(mean = 0, var = 0.1),
                       beta = list(mean = 0, var = 0.1))

beta_min <- -5
beta_max <- 5
beta_step <- 0.001

out <- bayesmr_noclus_optim(data, prior,
                            start = rnorm(2, mean = 0, sd = 20),
                            maxiter = 1000, tol = 1e-10,
                            beta_min = beta_min, beta_max = beta_min,
                            beta_step = beta_step,
                            n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                            eps_small = 1e-8, verbose = TRUE)
gamma_beta_post(out$gamma_best, out$beta_best, data, prior, log = TRUE, verbose = FALSE)

### beta full conditional

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

prior <- bayesmr_prior(gammaj = list(psi2 = 1),
                       Gammaj = list(tau2 = 1),
                       gamma = list(mean = 0, var = 1),
                       beta = list(mean = 0, var = 1))

beta_min <- -10
beta_max <- 10
beta_step <- 0.01
beta_seq <- seq(beta_min, beta_max, beta_step)
gamma_val <- 0.1
out <- logpost_beta(beta_seq, gamma, prior, data, log = FALSE)
plot(beta_seq, out, type = "l")
integrate(f = logpost_beta, lower = beta_min, upper = beta_max,
          subdivisions = 1e6, prior = prior, data = data, gamma = gamma_val,
          log = FALSE)

log_h2 <- Gj <- numeric(length(beta_seq))
for (i in 1:length(beta_seq)) {
  h2_j <- beta_seq[i]^2*prior[["gammaj"]][["psi2"]] +
            data$se_outcome^2 + prior[["Gammaj"]][["tau2"]]
  log_h2[i] <- -0.5*sum(log(h2_j))
  Gj[i] <- -0.5*sum((data$beta_outcome - beta_seq[i]*gamma_val)^2/h2_j)
}
plot(beta_seq, log_h2, type = "l")
lines(beta_seq, Gj, lty = 2)
abline(h = gamma_val^2/prior[["gammaj"]][["psi2"]])
abline(v = 0, lty = 2)
