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

prm.prop <- list(beta = 1.5)
burnin <- 20000
nsim <- 10000
seed <- 2301

set.seed(seed)

control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                random.start = TRUE, verbose = TRUE, nchains = 1, thin = 10,
                store.burnin = TRUE, threads = 1, parallel = "snow")
prior <- bayesmr_prior(gammaj = list(psi2 = mean(data[, 3]^2)),
                       Gammaj = list(tau2 = 1e-3),
                       gamma = list(mean = mean(data[, 1]),
                                    var = var(data[, 1])/nrow(data)),
                       beta = list(mean = 0.3, var = 0.2^2))
debugonce(bayesmr_fit)
sim_BMR <- bayesmr(zhaodata, control = control, prior = prior)
