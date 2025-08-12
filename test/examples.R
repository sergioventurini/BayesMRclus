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

# simulation setup
prm.prop <- list(beta = 1.5)
burnin <- 100000
nsim <- 200000

seed <- 2301
set.seed(seed)

control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                random.start = TRUE, verbose = TRUE, nchains = 3, thin = 100,
                store.burnin = TRUE, threads = 1, parallel = "snow")

prior <- bayesmr_prior(gammaj = list(psi2 = 0.1),
                       Gammaj = list(tau2 = 0.0001),
                       gamma = list(mean = 0, var = 0.1),
                       beta = list(mean = 0, var = 1.62))

# MCMC simulation
res_BMR <- bayesmr(zhaodata, control = control, prior = prior)

# MCMC summaries/graphs
summary(res_BMR)
plot(res_BMR, what = "trace", regex_pars = "beta")
plot(res_BMR, what = "trace", regex_pars = "gamma")
# plot(res_BMR, what = "hist_by_chain", regex_pars = "beta")
# plot(res_BMR, what = "dens_chains", regex_pars = "beta")
# plot(res_BMR, what = "intervals", regex_pars = "beta")
# plot(res_BMR, what = "acf_bar", regex_pars = "beta")

# ch <- 1
# beta_c <- res_BMR@results[[ch]]@beta.chain[1001:3000]
# plot(beta_c, type = "l")
# gamma_c <- res_BMR@results[[ch]]@gamma.chain[1001:3000]
# plot(gamma_c, type = "l")
# beta_gamma_c <- beta_c*gamma_c
# plot(beta_gamma_c, type = "l")
# summary(beta_gamma_c)
# plot(beta_c, gamma_c, pch = 20)

# coda analyses
library(coda)
res_BMR_sub <- subset(res_BMR, regex_pars = c("gamma", "beta"))
plot(res_BMR_sub)
cumuplot(res_BMR_sub)
gelman.plot(res_BMR_sub)
geweke.plot(res_BMR_sub)
raftery.diag(res_BMR_sub)
HPDinterval(res_BMR_sub)
heidel.diag(res_BMR_sub)
densplot(res_BMR_sub)
