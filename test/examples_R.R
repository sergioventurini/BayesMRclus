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

beta_z <- 0
beta_se_z <- 0.2
tau2_z <- 1e-5
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
tune <- .2
thin <- 10

# seed <- 2301
# set.seed(seed)

res_R <- mcmc_bayesmr(data_tmp, prior, iter, start, tune, proposal = "norm")
summary(res_R$draws[seq(burnin + 1, iter, by = thin), "gamma"])
plot(res_R$draws[seq(burnin + 1, iter, by = thin), "gamma"], type = "l",
     ylab = "gamma draws")
summary(res_R$draws[seq(burnin + 1, iter, by = thin), "beta"])
plot(res_R$draws[seq(burnin + 1, iter, by = thin), "beta"], type = "l",
     ylab = "beta draws")

### Numerical optimization of the joint posterior
map <- bayesmr_noclus_optim(data_tmp, prior, init = rnorm(2, mean = 0, sd = 20))
gamma_beta_post(map$par["gamma"], map$par["beta"],
                data_tmp, prior, log = TRUE)

### Wald ratios ###
data("bmi_sbp", package = "BayesMRclus")
data <- data.frame(beta_exposure = bmi_sbp[, "beta.exposure"],
                   beta_outcome = bmi_sbp[, "beta.outcome"],
                   se_exposure = bmi_sbp[, "se.exposure"],
                   se_outcome = bmi_sbp[, "se.outcome"])
betaj <- data_tmp$beta_outcome/data_tmp$beta_exposure
summary(betaj)
plot(betaj, pch = 19)
summary(data_tmp)
cm <- colMeans(data_tmp)
cm[2]/cm[1]  # these tend to be highly biased, especially
             # when the denominator is close to zero

library(MendelianRandomization)
mr_ivw(mr_input(data_tmp$beta_exposure, data_tmp$se_exposure,
                data_tmp$beta_outcome, data_tmp$se_outcome))
mr_median(mr_input(data_tmp$beta_exposure, data_tmp$se_exposure,
                   data_tmp$beta_outcome, data_tmp$se_outcome))

lm(beta_outcome ~ beta_exposure - 1, data = data_tmp, weights = se_outcome^-2)$coef
