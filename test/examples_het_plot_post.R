# (unnormalized) log-posterior contour plots showing the sampled values

library(patchwork)

### SETUP ###

# maximum a posteriori (MAP)
map <- bayesmr_noclus_het_optim(data_tmp, prior,
                                init = c(rnorm(2, mean = 0, sd = 5),
                                         runif(2, min = 0, max = 5)),
                                control = list(fnscale = -1, factr = 1e7,
                                               maxit = 1000, trace = 0))
map$par
map$message
# all.equal(map$value,
#           gamma_beta_psi_tau_post(map$par["gamma"], map$par["beta"],
#                                   map$par["psi"], map$par["tau"],
#                                   data_tmp, prior, log = TRUE))

gamma_min <- -.05
gamma_max <- 0.1
beta_min  <- -2.5
beta_max  <- 3
psi_min <- 0
psi_max <- 0.075
tau_min  <- 0
tau_max  <- 0.125

res <- 100
gamma_vals <- seq(gamma_min, gamma_max, length.out = res)
beta_vals  <- seq(beta_min, beta_max, length.out = res)
psi_vals <- seq(psi_min, psi_max, length.out = res)
tau_vals  <- seq(tau_min, tau_max, length.out = res)

beta_val <- map$par["beta"]
gamma_val <- map$par["gamma"]
psi_val <- map$par["psi"]
tau_val <- map$par["tau"]

res_BMR_sub <- subset(res_BMR, regex_pars = c("gamma", "beta", "psi", "tau"))
sample_points <- res_BMR_sub[[1]]
for (c in 2:control$nchains) {
  sample_points <- rbind(sample_points, res_BMR_sub[[c]])
}
sample_points <- data.frame(sample_points)
sample_points <- subset(
  sample_points,
  gamma >= gamma_min & gamma <= gamma_max &
    beta >= beta_min & beta <= beta_max &
    psi >= psi_min & psi <= psi_max &
    tau >= tau_min & tau <= tau_max
)

fill_scale <- scale_fill_viridis_d(
  option = "C",
  name = "Posterior"
)
axis_theme <- theme(
  axis.title.x = element_text(margin = margin(t = 6)),
  axis.title.y = element_text(margin = margin(r = 6)),
  axis.text.x  = element_text(),
  axis.text.y  = element_text()
)
margin_theme <- theme(
  plot.margin = margin(5.5, 5.5, 5.5, 5.5)
)

### PLOTS GENERATION ###

## gamma/beta
post_vals <- gamma_beta_psi_tau_post(gamma_vals, beta_vals,
                                     psi_val, tau_val,
                                     data_tmp, prior,
                                     log = TRUE, relative = TRUE)

df_plot <- expand.grid(gamma = gamma_vals, beta = beta_vals)
df_plot$posterior <- as.vector(post_vals)

p1 <- ggplot(df_plot, aes(x = gamma, y = beta, z = posterior)) +
  geom_contour_filled(bins = 20) +
  fill_scale +
  # geom_point(data = sample_points, aes(x = gamma, y = beta),
  #            color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = map$par[1], xend = map$par[1],
                   y = beta_min, yend = beta_max), color = "#21908C") +
  geom_segment(aes(x = gamma_min, xend = gamma_max,
                   y = map$par[2], yend = map$par[2]), color = "#21908C") +
  geom_point(x = map$par[1], y = map$par[2], color = "#21908C", size = 3) +
  labs(x = expression(gamma), y = expression(beta), fill = "Posterior") +
  coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(beta_min, beta_max)) +
  theme_minimal(base_size = 14) +
  axis_theme + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  margin_theme

## gamma/psi
post_vals <- gamma_beta_psi_tau_post(gamma_vals, beta_val,
                                     psi_vals, tau_val,
                                     data_tmp, prior,
                                     log = TRUE, relative = TRUE)

df_plot <- expand.grid(gamma = gamma_vals, psi = psi_vals)
df_plot$posterior <- as.vector(post_vals)

p2 <- ggplot(df_plot, aes(x = gamma, y = psi, z = posterior)) +
  geom_contour_filled(bins = 20) +
  fill_scale +
  # geom_point(data = sample_points, aes(x = gamma, y = psi),
  #            color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = map$par[1], xend = map$par[1],
                   y = psi_min, yend = psi_max), color = "#21908C") +
  geom_segment(aes(x = gamma_min, xend = gamma_max,
                   y = map$par[3], yend = map$par[3]), color = "#21908C") +
  geom_point(x = map$par[1], y = map$par[3], color = "#21908C", size = 3) +
  labs(x = expression(gamma), y = expression(psi), fill = "Posterior") +
  coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(psi_min, psi_max)) +
  theme_minimal(base_size = 14) +
  axis_theme + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  margin_theme

## gamma/tau
post_vals <- gamma_beta_psi_tau_post(gamma_vals, beta_val,
                                     psi_val, tau_vals,
                                     data_tmp, prior,
                                     log = TRUE, relative = TRUE)

df_plot <- expand.grid(gamma = gamma_vals, tau = tau_vals)
df_plot$posterior <- as.vector(post_vals)

p3 <- ggplot(df_plot, aes(x = gamma, y = tau, z = posterior)) +
  geom_contour_filled(bins = 20) +
  fill_scale +
  # geom_point(data = sample_points, aes(x = gamma, y = tau),
  #            color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = map$par[1], xend = map$par[1],
                   y = tau_min, yend = tau_max), color = "#21908C") +
  geom_segment(aes(x = gamma_min, xend = gamma_max,
                   y = map$par[4], yend = map$par[4]), color = "#21908C") +
  geom_point(x = map$par[1], y = map$par[4], color = "#21908C", size = 3) +
  labs(x = expression(gamma), y = expression(tau), fill = "Posterior") +
  coord_cartesian(xlim = c(gamma_min, gamma_max), ylim = c(tau_min, tau_max)) +
  theme_minimal(base_size = 14) +
  axis_theme + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  margin_theme

## beta/psi
post_vals <- gamma_beta_psi_tau_post(gamma_val, beta_vals,
                                     psi_vals, tau_val,
                                     data_tmp, prior,
                                     log = TRUE, relative = TRUE)

df_plot <- expand.grid(beta = beta_vals, psi = psi_vals)
df_plot$posterior <- as.vector(post_vals)

p4 <- ggplot(df_plot, aes(x = beta, y = psi, z = posterior)) +
  geom_contour_filled(bins = 20) +
  fill_scale +
  # geom_point(data = sample_points, aes(x = beta, y = psi),
  #            color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = map$par[2], xend = map$par[2],
                   y = psi_min, yend = psi_max), color = "#21908C") +
  geom_segment(aes(x = beta_min, xend = beta_max,
                   y = map$par[3], yend = map$par[3]), color = "#21908C") +
  geom_point(x = map$par[2], y = map$par[3], color = "#21908C", size = 3) +
  labs(x = expression(beta), y = expression(psi), fill = "Posterior") +
  coord_cartesian(xlim = c(beta_min, beta_max), ylim = c(psi_min, psi_max)) +
  theme_minimal(base_size = 14) +
  axis_theme + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  margin_theme

## beta/tau
post_vals <- gamma_beta_psi_tau_post(gamma_val, beta_vals,
                                     psi_val, tau_vals,
                                     data_tmp, prior,
                                     log = TRUE, relative = TRUE)

df_plot <- expand.grid(beta = beta_vals, tau = tau_vals)
df_plot$posterior <- as.vector(post_vals)

p5 <- ggplot(df_plot, aes(x = beta, y = tau, z = posterior)) +
  geom_contour_filled(bins = 20) +
  fill_scale +
  # geom_point(data = sample_points, aes(x = beta, y = tau),
  #            color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = map$par[2], xend = map$par[2],
                   y = tau_min, yend = tau_max), color = "#21908C") +
  geom_segment(aes(x = beta_min, xend = beta_max,
                   y = map$par[4], yend = map$par[4]), color = "#21908C") +
  geom_point(x = map$par[2], y = map$par[4], color = "#21908C", size = 3) +
  labs(x = expression(beta), y = expression(tau), fill = "Posterior") +
  coord_cartesian(xlim = c(beta_min, beta_max), ylim = c(tau_min, tau_max)) +
  theme_minimal(base_size = 14) +
  axis_theme + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  margin_theme

## psi/tau
post_vals <- gamma_beta_psi_tau_post(gamma_val, beta_val,
                                     psi_vals, tau_vals,
                                     data_tmp, prior,
                                     log = TRUE, relative = TRUE)

df_plot <- expand.grid(psi = psi_vals, tau = tau_vals)
df_plot$posterior <- as.vector(post_vals)

p6 <- ggplot(df_plot, aes(x = psi, y = tau, z = posterior)) +
  geom_contour_filled(bins = 20) +
  fill_scale +
  # geom_point(data = sample_points, aes(x = psi, y = tau),
  #            color = "gray", size = 0.1, inherit.aes = FALSE) +
  geom_segment(aes(x = map$par[3], xend = map$par[3],
                   y = tau_min, yend = tau_max), color = "#21908C") +
  geom_segment(aes(x = psi_min, xend = psi_max,
                   y = map$par[4], yend = map$par[4]), color = "#21908C") +
  geom_point(x = map$par[3], y = map$par[4], color = "#21908C", size = 3) +
  labs(x = expression(psi), y = expression(tau), fill = "Posterior") +
  coord_cartesian(xlim = c(psi_min, psi_max), ylim = c(tau_min, tau_max)) +
  theme_minimal(base_size = 14) +
  axis_theme + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  margin_theme

# # ---- Ridge curve implied by weak identifiability ----
# # constant C of the ridge tau^2 + beta^2 psi^2 = C
# C <- tau_val^2 + (beta_val^2) * psi_val^2
# psi_ridge <- seq(psi_min, psi_max, length.out = res)
# tau_ridge <- sqrt(pmax(0, C - (beta_val^2) * psi_ridge^2))
# ridge_df <- data.frame(psi = psi_ridge, tau = tau_ridge)
# 
# # add to plot
# p6 <- p6 +
#   geom_path(
#     data = ridge_df,
#     aes(x = psi, y = tau),
#     colour = "black",
#     linewidth = 1,
#     linetype = "longdash",
#     inherit.aes = FALSE
#   )
# p6

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

### PLOTS COMBINATION ###

# Add axis labels only to outer plots
# Top plot (p1): needs y-axis label
p1 <- p1 + theme(axis.title.y = element_text(margin = margin(r = 6)))

# Left column plots (p2, p3): need y-axis labels
p2 <- p2 + theme(axis.title.y = element_text(margin = margin(r = 6)))
p3 <- p3 + theme(axis.title.y = element_text(margin = margin(r = 6)))

# Bottom row plots (p3, p5, p6): need x-axis labels
p3 <- p3 + theme(axis.title.x = element_text(margin = margin(t = 6)))
p5 <- p5 + theme(axis.title.x = element_text(margin = margin(t = 6)))
p6 <- p6 + theme(axis.title.x = element_text(margin = margin(t = 6)))

layout <- "
A##
BD#
CEF
"

final_plot <-
  (p1 + p2 + p3 + p4 + p5 + p6) +
  plot_layout(
    design  = layout,
    guides  = "collect",
    widths  = c(1, 1, 1),
    heights = c(1, 1, 1)
  ) &
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.key.width = unit(2, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
final_plot

ggsave(
  "/Users/Sergio/Documents/Dati_VENTURINI/2_Research/1_Methods/BayesMRclus/Figs/fig1_lopost.png",
  final_plot,
  width = 15,
  height = 12,
  dpi = 300
)
