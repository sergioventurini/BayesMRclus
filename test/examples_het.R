library(BayesMRclus)

# prepare data
data("bmi_sbp", package = "BayesMRclus")
bmi_sbp <- subset(bmi_sbp, pval.selection < 5e-4)
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])

n <- nrow(data)
zhaodata <- new("bayesmr_data", data = data, n = n, harmonization = TRUE)
# summary(zhaodata)
# plot(zhaodata, se = TRUE)

# simulation setup
# note: psi and tau proposal standard deviations are on log scale (see paper)
prm.prop <- list(beta = .4, psi = .15, tau = .15)
burnin <- 100000
nsim <- 200000

# seed <- 2301
# set.seed(seed)

nchains <- 3
control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                psi.prop = prm.prop[["psi"]], tau.prop = prm.prop[["tau"]],
                random.start = TRUE, verbose = TRUE, nchains = nchains, thin = 50,
                store.burnin = TRUE, threads = ifelse(
                  nchains <= parallel::detectCores(),
                  nchains, parallel::detectCores() - 1),
                parallel = "snow")

alpha_psi <- 0.05 #alpha_eb(se = data_tmp[, "se_exposure"])
alpha_tau <- 2*alpha_psi #alpha_eb(se = data_tmp[, "se_outcome"])
prior <- bayesmr_prior(psi = list(alpha = alpha_psi, nu = 3),
                       tau = list(alpha = alpha_tau, nu = 3),
                       gamma = list(mean = 0, var = 1e2),
                       beta = list(mean = 0, var = 1e2))

# MCMC simulation
res_BMR <- bayesmr_het(zhaodata, control = control, prior = prior)

# MCMC summaries/graphs
summary(res_BMR)
plot(res_BMR, what = "trace", regex_pars = "psi")
plot(res_BMR, what = "trace", regex_pars = "tau")
plot(res_BMR, what = "trace", regex_pars = "beta")
plot(res_BMR, what = "trace", regex_pars = "gamma")
# plot(res_BMR, what = "hist_by_chain", regex_pars = "tau")
# plot(res_BMR, what = "dens_chains", regex_pars = "beta")
# plot(res_BMR, what = "intervals", regex_pars = "beta")
# plot(res_BMR, what = "acf_bar", regex_pars = "tau")

# coda analyses
library(coda)
# res_BMR_sub <- subset(res_BMR, regex_pars = c("gamma", "beta"))
res_BMR_sub <- subset(res_BMR, regex_pars = c("psi", "tau"))
plot(res_BMR_sub)
cumuplot(res_BMR_sub)
gelman.plot(res_BMR_sub)
geweke.plot(res_BMR_sub)
raftery.diag(res_BMR_sub)
HPDinterval(res_BMR_sub)
heidel.diag(res_BMR_sub)
densplot(res_BMR_sub)
autocorr.diag(res_BMR_sub)
effectiveSize(res_BMR_sub)
