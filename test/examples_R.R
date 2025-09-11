library(BayesMRclus)

# prepare data
data("bmi_sbp", package = "BayesMRclus")
bmi_sbp <- subset(bmi_sbp, pval.selection < 5e-4)
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

data_tmp <- data
flip <- data_tmp$beta_exposure < 0
data_tmp$beta_exposure[flip] <- -data_tmp$beta_exposure[flip]
data_tmp$beta_outcome[flip]  <- -data_tmp$beta_outcome[flip]

n <- nrow(data)
zhaodata <- new("bayesmr_data", data = data_tmp, n = n)
summary(zhaodata)
plot(zhaodata)

beta_z <- 0.3
beta_se_z <- 0.2
tau2_z <- 1e-3
prior <- bayesmr_prior(gammaj = list(psi2 = mean(data_tmp[, 3]^2)),
                       Gammaj = list(tau2 = tau2_z),
                       gamma = list(mean = mean(data_tmp[, 1]),
                                    var = var(data_tmp[, 1])/nrow(data_tmp)),
                       beta = list(mean = beta_z,
                                   var = beta_se_z^2))

# simulation setup
burnin <- 10000
nsim <- 20000

iter <- burnin + nsim
start <- rnorm(2)
tune <- .65

seed <- 2301
set.seed(seed)

res_R <- mcmc_bayesmr(data_tmp, prior, iter, start, tune, proposal = "norm")
summary(res_R$draws[, 1])
plot(res_R$draws[, 1], type = "l")
summary(res_R$draws[, 2])
plot(res_R$draws[, 2], type = "l")

### Numerical optimization of the joint posterior

beta_min <- -5
beta_max <- 5
beta_step <- 0.001

out <- bayesmr_noclus_optim(data_tmp, prior,
                            start = rnorm(2, mean = 0, sd = 20),
                            maxiter = 1000, tol = 1e-10,
                            beta_min = beta_min, beta_max = beta_min,
                            beta_step = beta_step,
                            n = 1000, tol_x = 1e-10, tol_f = 1e-12,
                            eps_small = 1e-8, verbose = TRUE)
gamma_beta_post(out$gamma_best, out$beta_best, data, prior, log = TRUE, verbose = FALSE)

### beta full conditional

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
          subdivisions = 1e6, prior = prior, data = data_tmp, gamma = gamma_val,
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

### Wald ratios ###

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])
betaj <- data$beta_outcome/data$beta_exposure
summary(betaj)
plot(betaj, pch = 19)
summary(data)
cm <- colMeans(data)
cm[2]/cm[1]

library(MendelianRandomization)
mr_ivw(mr_input(data$beta_exposure, data$se_exposure,
                data$beta_outcome, data$se_outcome))
mr_median(mr_input(data$beta_exposure, data$se_exposure,
                   data$beta_outcome, data$se_outcome))

lm(beta_outcome ~ beta_exposure - 1, data = data, weights = se_outcome^-2)$coef

### profile log-likelihood

library(BayesMRclus)

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

psi2 <- .7
tau2 <- .001
beta_min <- -10
beta_max <- 10
beta_step <- 0.01
beta_seq <- seq(beta_min, beta_max, beta_step)
llik <- profile_loglik_all(beta_seq, psi2, tau2, data)
plot(beta_seq, llik, type = "l")
abline(v = 0.029)

beta <- .03
tau2 <- .0001
psi2_min <- 0
psi2_max <- 5
psi2_step <- 0.001
psi2_seq <- seq(psi2_min, psi2_max, psi2_step)
llik <- profile_loglik_all(beta, psi2_seq, tau2, data)
plot(psi2_seq, llik, type = "l")
abline(v = 0.7)

beta <- .03
psi2 <- .7
tau2_min <- 0
tau2_max <- 1
tau2_step <- 0.001
tau2_seq <- seq(tau2_min, tau2_max, tau2_step)
llik <- profile_loglik_all(beta, psi2, tau2_seq, data)
plot(tau2_seq, llik, type = "l")
abline(v = 0)

out <- fit_profile_loglik(data)
out
cov_par <- solve(out$opt$hessian[1:2, 1:2])
sqrt(cov_par[1, 1])
out$psi2 * sqrt(cov_par[2, 2])

# another data set
data("amd_hdl", package = "BayesMRclus")
data <- data.frame(beta_exposure = amd_hdl[, "beta_XL.HDL.C"],
                   beta_outcome = amd_hdl[, "beta_amd"],
                   se_exposure = amd_hdl[, "se_XL.HDL.C"],
                   se_outcome = amd_hdl[, "se_amd"])

psi2 <- .57
tau2 <- .001
beta_min <- -2
beta_max <- 2
beta_step <- 0.01
beta_seq <- seq(beta_min, beta_max, beta_step)
llik <- profile_loglik_all(beta_seq, psi2, tau2, data)
plot(beta_seq, llik, type = "l")
abline(v = 0.045)

beta <- .045
tau2 <- .0001
psi2_min <- 0
psi2_max <- 3
psi2_step <- 0.001
psi2_seq <- seq(psi2_min, psi2_max, psi2_step)
llik <- profile_loglik_all(beta, psi2_seq, tau2, data)
plot(psi2_seq, llik, type = "l")
abline(v = 0.57)

beta <- .045
psi2 <- .57
tau2_min <- 0
tau2_max <- 1
tau2_step <- 0.001
tau2_seq <- seq(tau2_min, tau2_max, tau2_step)
llik <- profile_loglik_all(beta, psi2, tau2_seq, data)
plot(tau2_seq, llik, type = "l")
abline(v = 0)

out <- fit_profile_loglik(data)
out

data_tmp <- data
flip <- data_tmp$beta_exposure < 0
data_tmp$beta_exposure[flip] <- -data_tmp$beta_exposure[flip]
data_tmp$beta_outcome[flip]  <- -data_tmp$beta_outcome[flip]
plot(data_tmp[, 1:2], pch = 19, ylim = c(-.2, .2))
# summary(lm(beta_outcome ~ beta_exposure - 1, data = data_tmp))

library(mr.raps)
fit_ps <- mr.raps.simple(b_out  = data_tmp$beta_outcome,
                         se_out = data_tmp$se_outcome,
                         b_exp  = data_tmp$beta_exposure,
                         se_exp = data_tmp$se_exposure)

fit_raps <- mr.raps.overdispersed.robust(b_out  = data_tmp$beta_outcome,
                                         se_out = data_tmp$se_outcome,
                                         b_exp  = data_tmp$beta_exposure,
                                         se_exp = data_tmp$se_exposure)

plot(data_tmp$beta_exposure, data_tmp$beta_outcome, pch = 19,
     xlab = "SNP effect on HDL", ylab = "SNP effect on AMD")
with(data_tmp, {
  # horizontal (exposure) error bars
  segments(beta_exposure - se_exposure, beta_outcome,
           beta_exposure + se_exposure, beta_outcome)
  # vertical (outcome) error bars
  segments(beta_exposure, beta_outcome - se_outcome,
           beta_exposure, beta_outcome + se_outcome)
})
abline(a = 0, b = fit_ps$beta.hat,   lwd = 2)           # "Simple" PS
abline(a = 0, b = fit_raps$beta.hat, lwd = 2, lty = 2)  # Overdispersion + robust
legend("bottomleft",
       c("Simple", "Overdispersion + Robust loss"),
       lwd = 2, lty = c(1,2))

### MCMC using profile likelihood

data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

data_tmp <- data
# flip <- data_tmp$beta_exposure < 0
# data_tmp$beta_exposure[flip] <- -data_tmp$beta_exposure[flip]
# data_tmp$beta_outcome[flip]  <- -data_tmp$beta_outcome[flip]
# plot(data_tmp[, 1:2], pch = 19, ylim = c(-.2, .2))

beta_z <- 0.3
beta_se_z <- 0.2
prec_psi2_shape <- 1e-3
prec_psi2_rate <- 1e-3
prec_tau2_shape <- 1e-3
prec_tau2_rate <- 1e-3
prior <- list(beta = list(mean = beta_z,
                          var = beta_se_z^2),
              prec_psi2 = list(shape = prec_psi2_shape,
                               rate = prec_psi2_rate),
              prec_tau2 = list(shape = prec_tau2_shape,
                               rate = prec_tau2_rate)
              )

# simulation setup
burnin <- 10000
nsim <- 20000

iter <- burnin + nsim
start <- c(0.5, 0.3, 0.0001) #c(rnorm(1), runif(2, 0.0001, 0.1))
tune <- c(.15, 100, 10000)

seed <- 2301
set.seed(seed)

res_R <- mcmc_bayesmr_profile(data_tmp, prior, iter, start, tune, verbose = TRUE)
summary(res_R$draws[, 1])
plot(res_R$draws[, 1], type = "l")
summary(res_R$draws[, 2])
plot(res_R$draws[, 2], type = "l")
summary(res_R$draws[, 3])
plot(res_R$draws[, 3], type = "l")
res_R$accept_rate
