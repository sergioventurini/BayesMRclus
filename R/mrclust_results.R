#' Run MR-Clust and summarise results.
#'
#' Fits the MR-Clust EM algorithm ([mrclust::mr_clust_em()]) to SNP-level
#' summary data and returns two diagnostic plots and a cluster summary table.
#'
#' @param data A data frame with columns `beta.exposure`, `beta.outcome`,
#'   `se.exposure`, `se.outcome`, and `SNP`.
#' @param title A character string used in the plot titles. Default is
#'   `"X -> Y"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{`p1`}{A `ggplot2` weighted-density plot coloured by cluster
#'     assignment.}
#'   \item{`p2`}{A `ggplot2` scatter plot of Wald ratio versus precision,
#'     coloured by cluster.}
#'   \item{`tab`}{A data frame summarising each cluster: cluster mean, number of
#'     SNPs, mean assignment probability, and number of high-confidence
#'     assignments.}
#'   \item{`cluster`}{A data frame with columns `SNP`, `cluster`, and
#'     `cluster_class` giving per-SNP cluster memberships.}
#' }
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso [mrclust::mr_clust_em()].
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' res <- mrclust_results(ldl_cad, title = "LDL-C -> CAD")
#' print(res$p1)
#' print(res$p2)
#' print(res$tab)
#' }
#'
#' @importFrom mrclust mr_clust_em
#' @export
mrclust_results <- function(data, title = "X -> Y") {
  wald    <- data$beta.outcome / data$beta.exposure
  wald_se <- data$se.outcome   / abs(data$beta.exposure)

  suppressWarnings(
    res <- mrclust::mr_clust_em(theta = wald, theta_se = wald_se,
                                bx = data$beta.exposure, bxse = data$se.exposure,
                                by = data$beta.outcome,  byse = data$se.outcome,
                                obs_names = data$SNP)
  )

  best <- res$results$best

  # Cluster mean lines, one per cluster_class.
  clust_lines <- dplyr::select(
    dplyr::slice(
      dplyr::group_by(best, .data$cluster_class),
      1
    ),
    .data$cluster_class,
    .data$cluster_mean
  )

  # Plot 1: weighted density by cluster.
  p1 <- ggplot2::ggplot(best, ggplot2::aes(x = theta, fill = cluster_class,
                        weight = 1 / theta_se^2)) +
    ggplot2::geom_density(alpha = 0.45, bw = "SJ") +
    ggplot2::geom_vline(data = clust_lines,
                        ggplot2::aes(xintercept = cluster_mean, colour = cluster_class),
                        linetype = "dashed", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    ggplot2::labs(x     = "Wald ratio (SNP-specific causal estimate)",
                  y     = "Weighted density",
                  fill  = "Cluster",
                  colour = "Cluster mean",
                  title = paste0("MR-Clust: ", title, " cluster structure")) +
    ggplot2::theme_minimal()

  # Plot 2: scatter plot coloured by cluster and sized by precision.
  p2 <- ggplot2::ggplot(best, ggplot2::aes(x = theta, y = 1 / theta_se,
                        colour = cluster_class,
                        linewidth = 1 / theta_se^2)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = theta - 1.96 * theta_se,
                                         xmax = theta + 1.96 * theta_se),
                                         width = 0, alpha = 0.3, linewidth = 0.4) +
    ggplot2::geom_vline(data = clust_lines,
                        ggplot2::aes(xintercept = cluster_mean, colour = cluster_class),
                        linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    ggplot2::scale_size_continuous(guide = "none") +
    ggplot2::labs(x      = "Wald ratio",
                  y      = "Precision (1 / SE)",
                  colour = "Cluster",
                  title  = paste0("MR-Clust: SNP assignments (", title, ")")) +
    ggplot2::theme_minimal()

  # Cluster summary table.
  tab <- dplyr::arrange(
    dplyr::summarise(
      dplyr::group_by(best, .data$cluster_class, .data$cluster_mean),
      n_snps = dplyr::n(),
      mean_prob = round(mean(.data$probability), 3),
      n_high_conf = sum(.data$probability > 0.8),
      .groups = "drop"
    ),
    .data$cluster_mean
  )

  list(
    p1 = p1,
    p2 = p2,
    tab = tab,
    cluster = data.frame(
      SNP = rownames(data),
      cluster = best$cluster,
      cluster_class = best$cluster_class
    )
  )
}
