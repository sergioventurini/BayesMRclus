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
tune <- 1.5

seed <- 2301
set.seed(seed)

res_R <- mcmc_bayesmr(data, hpar, iter, start, tune)
summary(res_R$draws[, 1])
plot(res_R$draws[, 1], type = "l")
summary(res_R$draws[, 2])
plot(res_R$draws[, 2], type = "l")
