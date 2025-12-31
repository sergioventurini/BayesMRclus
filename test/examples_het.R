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
prm.prop <- list(beta = .2, psi = .04, tau = .04)
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

# maximum a posteriori (MAP)
map <- bayesmr_noclus_het_optim(data_tmp, prior,
                                init = c(rnorm(2, mean = 0, sd = 5),
                                         runif(2, min = 0, max = 5)),
                                control = list(fnscale = -1, factr = 1e7,
                                               maxit = 1000, trace = 0))
map$par
map$message
all.equal(map$value,
          gamma_beta_psi_tau_post(map$par["gamma"], map$par["beta"],
                                  map$par["psi"], map$par["tau"],
                                  data_tmp, prior, log = TRUE))

# (unnormalized) joint psi/tau posterior contour plot showing the sampled values
psi_min <- 0
psi_max <- 0.035
tau_min  <- 0
tau_max  <- 0.075

res <- 100
psi_vals <- seq(psi_min, psi_max, length.out = res)
tau_vals  <- seq(tau_min, tau_max, length.out = res)

beta_val <- map$par["beta"]
gamma_val <- map$par["gamma"]

post_vals <- gamma_beta_psi_tau_post(gamma_val, beta_val,
                                     psi_vals, tau_vals,
                                     data_tmp, prior, log = TRUE)

df_plot <- expand.grid(psi = psi_vals, tau = tau_vals)
df_plot$posterior <- as.vector(post_vals)

res_BMR_sub <- subset(res_BMR, regex_pars = c("psi", "tau"))
sample_points <- res_BMR_sub[[1]]
if (control$nchains > 1) {
  for (c in 2:control$nchains) {
    sample_points <- rbind(sample_points, res_BMR_sub[[c]])
  }
}
sample_points <- data.frame(sample_points)
sample_points <- subset(
  sample_points,
  psi >= psi_min & psi <= psi_max &
    tau >= tau_min & tau <= tau_max
)

p <-
  ggplot(df_plot, aes(x = psi, y = tau, z = posterior)) +
    geom_contour_filled(bins = 20) +
    scale_fill_viridis_d(option = "C") +
    geom_point(data = sample_points, aes(x = psi, y = tau),
               color = "gray", size = 0.1, inherit.aes = FALSE) +
    geom_segment(aes(x = map$par[3], xend = map$par[3],
                     y = tau_min, yend = tau_max), color = "#21908C") +
    geom_segment(aes(x = psi_min, xend = psi_max,
                     y = map$par[4], yend = map$par[4]), color = "#21908C") +
    geom_point(x = map$par[3], y = map$par[4], color = "#21908C", size = 3) +
    labs(x = expression(psi), y = expression(tau), fill = "Posterior") +
    coord_cartesian(xlim = c(psi_min, psi_max), ylim = c(tau_min, tau_max)) +
    theme_minimal(base_size = 14)
p

# # ---- Ridge curve implied by weak identifiability ----
# # Posterior mode
# psi_star <- map$par["psi"]
# tau_star <- map$par["tau"]
# beta_est <- beta_val
# 
# # constant C of the ridge tau^2 + beta^2 psi^2 = C
# C <- tau_star^2 + (beta_est^2) * psi_star^2
# psi_ridge <- seq(psi_min, psi_max, length.out = 400)
# tau_ridge <- sqrt(pmax(0, C - (beta_est^2) * psi_ridge^2))
# ridge_df <- data.frame(psi = psi_ridge, tau = tau_ridge)
# 
# # add to plot
# p <- p +
#   geom_path(
#     data = ridge_df,
#     aes(x = psi, y = tau),
#     colour = "black",
#     linewidth = 1,
#     linetype = "longdash",
#     inherit.aes = FALSE
#   )
# p

# # interactive 3D plot
# library(plotly)
# zmat <- matrix(df_plot$posterior, nrow = length(psi_vals), ncol = length(tau_vals))
# plot_ly(
#   x = psi_vals,
#   y = tau_vals,
#   z = zmat
# ) |>
#   add_surface(colorscale = "Viridis") |>
#   layout(
#     scene = list(
#       xaxis = list(title = "psi"),
#       yaxis = list(title = "tau"),
#       zaxis = list(title = "log posterior")
#     )
#   )
