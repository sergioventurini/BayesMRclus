library(ggplot2)
library(mr.raps)
library(MRPATH)
data(hdl_chd)

plot_data <- subset(hdl_chd, pval.selection <= 5e-8)
sel <- plot_data$beta.exposure < 0
plot_data$beta.exposure[sel] <- -plot_data$beta.exposure[sel]
plot_data$beta.outcome[sel] <- -plot_data$beta.outcome[sel]

# Fit MR Egger regression
# Fit MR-Egger
weights <- 1 / plot_data$se.outcome^2
egger_fit <- lm(beta.outcome ~ beta.exposure, 
                data = plot_data, 
                weights = weights)
egger_slope <- egger_fit[[1]][2]
egger_intercept <- egger_fit[[1]][1]

# Fit MR Raps
# If not available, we'll use a weighted regression as approximation
raps_fit <- mr.raps(data = plot_data,
                    diagnostics = FALSE,
                    over.dispersion = FALSE,
                    loss.function = "huber")
raps_slope <- raps_fit$beta.hat
raps_intercept <- 0  # MR Raps assumes no pleiotropy

# Create the plot
p <- ggplot(plot_data, aes(x = beta.exposure, y = beta.outcome)) +
  # Add points with error bars
  geom_hline(yintercept = 0, color = "gray", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "gray", linewidth = 0.5) +
  geom_point(size = 1, alpha = 0.6) +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome, 
                    ymax = beta.outcome + se.outcome),
                alpha = 0.3, width = 0) +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure,
                     xmax = beta.exposure + se.exposure),
                 alpha = 0.3, height = 0) +
  # Add MR Egger line
  geom_abline(intercept = egger_intercept, slope = egger_slope,
              color = "#F8766D", linewidth = 0.8) +
  # Add MR Raps line
  geom_abline(intercept = raps_intercept, slope = raps_slope,
              color = "#00BFC4", linewidth = 0.8) +
  # Add confidence interval for MR Raps (approximate)
  # Labels and theme
  labs(x = "SNP association with HDL-C",
       y = "SNP association with CHD") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = c(0.85, 0.15),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "white", color = NA)
  ) +
  # Add legend manually
  annotate("text", x = Inf, y = -Inf, 
           label = "— MR Egger  — MR Raps",
           hjust = 1.05, vjust = -0.5, size = 3.5,
           color = "black")

# Display the plot
print(p)

# Optional: Save the plot
# ggsave("figure1a_reproduction.pdf", p, width = 6, height = 5)
