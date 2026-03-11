#' @importFrom mrclust mr_clust_em
#' @export
mrclust_results <- function(data, title = "X → Y") {
  wald    <- data$beta.outcome / data$beta.exposure
  wald_se <- data$se.outcome   / abs(data$beta.exposure)
  
  suppressWarnings(
    res <- mrclust::mr_clust_em(theta = wald, theta_se = wald_se,
                                bx = data$beta.exposure, bxse = data$se.exposure,
                                by = data$beta.outcome,  byse = data$se.outcome,
                                obs_names = data$SNP)
  )
  
  best <- res$results$best
  
  # Cluster mean lines (one per cluster_class)
  clust_lines <- best %>%
    group_by(cluster_class) %>%
    slice(1) %>%
    select(cluster_class, cluster_mean)
  
  # ── Plot 1: weighted density by cluster ───────────────────────────────────────
  p1 <- ggplot(best, aes(x = theta, fill = cluster_class,
                   weight = 1 / theta_se^2)) +
    geom_density(alpha = 0.45, bw = "SJ") +
    geom_vline(data = clust_lines,
               aes(xintercept = cluster_mean, colour = cluster_class),
               linetype = "dashed", linewidth = 0.9) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    labs(x     = "Wald ratio (SNP-specific causal estimate)",
         y     = "Weighted density",
         fill  = "Cluster",
         colour = "Cluster mean",
         title = paste0("MR-Clust: ", title, " cluster structure")) +
    theme_minimal()
  
  # ── Plot 2: scatter plot coloured by cluster, sized by precision ──────────────
  p2 <- ggplot(best, aes(x = theta, y = 1 / theta_se,
                   colour = cluster_class,
                   linewidth = 1 / theta_se^2)) +
    geom_point(alpha = 0.6) +
    geom_errorbarh(aes(xmin = theta - 1.96 * theta_se,
                       xmax = theta + 1.96 * theta_se),
                   width = 0, alpha = 0.3, linewidth = 0.4) +
    geom_vline(data = clust_lines,
               aes(xintercept = cluster_mean, colour = cluster_class),
               linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    scale_size_continuous(guide = "none") +
    labs(x      = "Wald ratio",
         y      = "Precision (1 / SE)",
         colour = "Cluster",
         title  = paste0("MR-Clust: SNP assignments (", title, ")")) +
    theme_minimal()
  
  # ── Cluster summary table ─────────────────────────────────────────────────────
  tab <- best %>%
    group_by(cluster_class, cluster_mean) %>%
    summarise(n_snps      = n(),
              mean_prob   = round(mean(probability), 3),
              n_high_conf = sum(probability > 0.8),
              .groups     = "drop") %>%
    arrange(cluster_mean)
  
  list(p1 = p1, p2 = p2, tab = tab)
}
