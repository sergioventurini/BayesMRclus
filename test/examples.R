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
# summary(zhaodata)
# plot(zhaodata)

# simulation setup
prm.prop <- list(beta = .1)
burnin <- 100000
nsim <- 200000

# seed <- 2301
# set.seed(seed)

nchains <- 3
control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                random.start = TRUE, verbose = TRUE, nchains = nchains, thin = 50,
                store.burnin = TRUE, threads = ifelse(
                  nchains <= parallel::detectCores(),
                  nchains, parallel::detectCores() - 1),
                parallel = "snow")

# tau2_est <- tau2_dl(data_tmp, secondorder = TRUE)
prior <- bayesmr_prior(gammaj = list(psi2 = .0001080392),
                       Gammaj = list(tau2 = .0005981383),
                       gamma = list(mean = 0, var = 1e2),
                       beta = list(mean = 0, var = 1e2))

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

# (unnormalized) posterior contour plot showing the sampled values
gamma_min <- 0
gamma_max <- 0.03
beta_min  <- -0.5
beta_max  <- 1

res <- 100
gamma_vals <- seq(gamma_min, gamma_max, length.out = res)
beta_vals  <- seq(beta_min, beta_max, length.out = res)

post_vals <- gamma_beta_post(gamma_vals, beta_vals,
                             data_tmp, prior, log = TRUE)

df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
df_plot$posterior <- as.vector(post_vals)

map <- bayesmr_noclus_optim(data_tmp, prior, init = rnorm(2, mean = 0, sd = 20))
all.equal(map$value,
          gamma_beta_post(map$par["gamma"], map$par["beta"],
                          data_tmp, prior, log = TRUE))

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
  geom_segment(aes(x = map$par[1], xend = map$par[1],
        y = beta_min, yend = beta_max), color = "#21908C") +
  geom_segment(aes(x = gamma_min, xend = gamma_max,
        y = map$par[2], yend = map$par[2]), color = "#21908C") +
  geom_point(x = map$par[1], y = map$par[2], color = "#21908C", size = 3) +
  labs(x = expression(gamma), y = expression(beta), fill = "Posterior") +
  coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
  theme_minimal(base_size = 14)

# log-likelihood contour plot
gamma_min <- 0
gamma_max <- 0.03
beta_min  <- -0.5
beta_max  <- 1

res <- 100
gamma_vals <- seq(gamma_min, gamma_max, length.out = res)
beta_vals  <- seq(beta_min, beta_max, length.out = res)

ll_vals <- bayesmr_logLik(gamma_vals, beta_vals, data_tmp, prior,
                          log = TRUE)

df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
df_plot$ll <- as.vector(ll_vals)

mle <- bayesmr_noclus_mle(data_tmp, prior, init = rnorm(2, mean = 0, sd = 20))
all.equal(mle$value,
          bayesmr_logLik(mle$par["gamma"], mle$par["beta"],
                         data_tmp, prior, log = TRUE))

ggplot(df_plot, aes(x = gamma, y = beta, z = ll)) +
  geom_contour_filled(bins = 20) +
  scale_fill_viridis_d(option = "C") +
  geom_segment(aes(x = mle$par[1], xend = mle$par[1],
                   y = beta_min, yend = beta_max), color = "#21908C") +
  geom_segment(aes(x = gamma_min, xend = gamma_max,
                   y = mle$par[2], yend = mle$par[2]), color = "#21908C") +
  geom_point(x = mle$par[1], y = mle$par[2], color = "#21908C", size = 3) +
  labs(x = expression(gamma), y = expression(beta), fill = "Log-likelihood") +
  coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
  theme_minimal(base_size = 14)
