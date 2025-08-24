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

prior <- bayesmr_prior(gammaj = list(psi2 = 0.01),
                       Gammaj = list(tau2 = 0.5),
                       gamma = list(mean = 0, var = 0.01),
                       beta = list(mean = 0, var = 1))

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

# contour plot showing the sampled values
gamma_min <- -0.05
gamma_max <- 0.05
beta_min  <- -4
beta_max  <- 4

res <- 100
gamma_vals <- seq(gamma_min, gamma_max, length.out = res)
beta_vals  <- seq(beta_min, beta_max, length.out = res)

post_vals <- gamma_beta_post(gamma_vals, beta_vals, data, prior, log = FALSE, verbose = FALSE)

df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
df_plot$posterior <- as.vector(post_vals)

glob_max <- bayesmr_noclus_optim(data, prior, start = rep(0, 2), maxiter = 1000, tol = 1e-10,
                                 beta_min = beta_min, beta_max = beta_max, beta_step = 0.001,
                                 n = 1000, tol_x = 1e-10, tol_f = 1e-12, eps_small = 1e-8,
                                 verbose = FALSE)

res_BMR_sub <- subset(res_BMR, regex_pars = c("gamma", "beta"))
sample_points <- res_BMR_sub[[1]]
for (c in 2:control$nchains) {
  sample_points <- rbind(sample_points, res_BMR_sub[[c]])
}
sample_points <- data.frame(sample_points)
sample_points <- subset(
  sample_points,
  gamma >= gamma_min & gamma <= gamma_max &
    beta >= beta_min & beta <= beta_max
)

ggplot(df_plot, aes(x = gamma, y = beta, z = posterior)) +
  geom_contour_filled(bins = 20) +
  scale_fill_viridis_d(option = "C") +
  geom_point(data = sample_points, aes(x = gamma, y = beta),
             color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = glob_max$gamma_best, xend = glob_max$gamma_best,
        y = beta_min, yend = beta_max), color = "#21908C") +
  geom_segment(aes(x = gamma_min, xend = gamma_max,
        y = glob_max$beta_best, yend = glob_max$beta_best), color = "#21908C") +
  geom_point(x = glob_max$gamma_best, y = glob_max$beta_best, color = "#21908C", size = 3) +
  labs(x = expression(gamma), y = expression(beta), fill = "Posterior") +
  coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
  theme_minimal(base_size = 14)

###

# weight of beta prior over beta full conditional
b <- seq(-5, 5, .01)
par <- list(gamma = 0)
prior <- bayesmr_prior(gammaj = list(psi2 = 0.01),
                       Gammaj = list(tau2 = 0.5),
                       gamma = list(mean = 0, var = 0.01),
                       beta = list(mean = 0, var = 1))
hpar <- list(mu_beta = prior$beta$mean,
             sigma2_beta = prior$beta$var,
             psi2 = prior$gammaj$psi2,
             tau2 = prior$Gammaj$tau2)

out <- logpost_beta_util(b, par, hpar, data)
beta_prior_w <- out[, 2]/rowSums(out)
beta_prior_w <- out[, 2] - out[, 1]
beta_prior_w <- out[, 2]/out[, 1]
plot(b, beta_prior_w, type = "l")

summary(out)
plot(b, out[, 1], type = "l", ylim = c(min(out), max(out)))
lines(b, out[, 2], col = "red")

summary(exp(out))
plot(b, exp(out[, 1]), type = "l", ylim = c(min(exp(out)), max(exp(out))))
lines(b, exp(out[, 2]), col = "red")
