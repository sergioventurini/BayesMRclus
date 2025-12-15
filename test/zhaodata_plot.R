library(BayesMRclus)

# prepare data
data("bmi_sbp", package = "BayesMRclus")
bmi_sbp <- subset(bmi_sbp, pval.selection < 5e-8)
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

library(ggplot2)

# Create 95% confidence intervals for error bars
zhaodata@data$lower_x <- zhaodata@data$beta_exposure - zhaodata@data$se_exposure
zhaodata@data$upper_x <- zhaodata@data$beta_exposure + zhaodata@data$se_exposure
zhaodata@data$lower_y <- zhaodata@data$beta_outcome - zhaodata@data$se_outcome
zhaodata@data$upper_y <- zhaodata@data$beta_outcome + zhaodata@data$se_outcome

# Two fitted models:
library(mr.raps)

zhaodata2 <- zhaodata
colnames(zhaodata2@data) <- c("beta.exposure", "beta.outcome",
                              "se.exposure", "se.outcome")

fit_simple <- mr.raps(zhaodata2@data, diagnostics = FALSE,
                      over.dispersion = FALSE, loss = "l2")

# Overdispersion + robust fit (Huber loss)
# Requires MASS for rlm()
fit_robust <- mr.raps(zhaodata2@data, diagnostics = FALSE,
                      over.dispersion = TRUE, loss = "huber")

# New dataframe for prediction lines
newdat <- data.frame(
  beta_exposure = seq(-.01,
                      .1,
                      length.out = 200)
)

# Predictions and 95% CI for both models
newdat$simple_fit  <- fit_simple$beta.hat*newdat$beta_exposure
newdat$simple_lwr  <- (fit_simple$beta.hat - 1.96*fit_simple$beta.se)*newdat$beta_exposure
newdat$simple_upr  <- (fit_simple$beta.hat + 1.96*fit_simple$beta.se)*newdat$beta_exposure

# For rlm, manual bootstrap or normal-approx CI; here we use normal approx:
newdat$robust_fit <- fit_robust$beta.hat*newdat$beta_exposure
newdat$robust_lwr <- (fit_robust$beta.hat - 1.96*fit_robust$beta.se)*newdat$beta_exposure
newdat$robust_upr <- (fit_robust$beta.hat + 1.96*fit_robust$beta.se)*newdat$beta_exposure

# Plot
ggplot(zhaodata@data, aes(x = beta_exposure, y = beta_outcome)) +
  
  # Error bars
  geom_errorbar(aes(ymin = lower_y, ymax = upper_y),
                width = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = lower_x, xmax = upper_x),
                 height = 0, colour = "grey50") +
  
  # Points
  geom_point(size = 2, shape = 21, fill = "black") +
  
  # Simple regression line + CI
  geom_line(data = newdat, aes(x = beta_exposure, y = simple_fit,
                               colour = "Simple")) +
  geom_line(data = newdat, aes(x = beta_exposure, y = simple_lwr,
                               colour = "Simple"), linetype = "dashed") +
  geom_line(data = newdat, aes(x = beta_exposure, y = simple_upr,
                               colour = "Simple"), linetype = "dashed") +
  
  # Robust regression line + CI
  geom_line(data = newdat, aes(x = beta_exposure, y = robust_fit,
                               colour = "Overdispersion + Robust loss")) +
  geom_line(data = newdat, aes(x = beta_exposure, y = robust_lwr,
                               colour = "Overdispersion + Robust loss"),
            linetype = "dashed") +
  geom_line(data = newdat, aes(x = beta_exposure, y = robust_upr,
                               colour = "Overdispersion + Robust loss"),
            linetype = "dashed") +
  
  scale_colour_manual(values = c("Simple" = "#C44E52",
                                 "Overdispersion + Robust loss" = "#4C9FAD")) +
  
  labs(x = "SNP effect on BMI",
       y = "SNP effect on SBP",
       colour = "Method") +
  
  scale_x_continuous(limits = c(-0.002, 0.095),
                       breaks = seq(0, 0.1, 0.025)) +
  scale_y_continuous(limits = c(-0.055, 0.0875),
                     breaks = seq(-0.05, 0.08, 0.05)) +

  theme_minimal(base_size = 16) +
  
  theme(legend.position = "bottom") +
  
  theme(panel.border = element_rect(colour = "black", fill = NA, size = .5))
