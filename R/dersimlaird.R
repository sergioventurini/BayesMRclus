#' Compute DerSimonian-Laird Tau-Squared Estimate
#'
#' Calculates the DerSimonian-Laird estimate of the between-study variance
#' (tau-squared) for meta-analysis of Mendelian randomization studies.
#'
#' @param data A matrix or data frame with 4 columns:
#'   \describe{
#'     \item{Column 1}{Estimated genetic associations with the exposure (\code{gammaj_hat})}
#'     \item{Column 2}{Estimated genetic associations with the outcome (\code{Gammaj_hat})}
#'     \item{Column 3}{Standard errors of exposure associations (\code{sigmaj_X})}
#'     \item{Column 4}{Standard errors of outcome associations (\code{sigmaj_Y})}
#'   }
#'   Each row represents a different genetic variant.
#' @param secondorder Length-one logical vector indicating whether to apply a second
#'   order delta method expansion for the variance of a ratio. Default to \code{TRUE}.
#' @param beta_max Length-one numeric vector representing the maximum Wald ratio value
#'   to use. That is, all the genetic variants whose Wald ratio will be (in absolute
#'   terms) larger than \code{beta_max} will be removed before computing the
#'   DerSimonian-Laird estimate.
#' @param se_max Length-one numeric vector representing the maximum Wald ratio standard
#'   error to use. That is, all the genetic variants whose Wald ratio standard error
#'   will be (in absolute terms) larger than \code{se_max} will be removed before
#'   computing the DerSimonian-Laird estimate.
#'
#' @return A list that contains the DerSimonian-Laird estimate of tau-squared, the Wald
#'   ratio estimates and the corresponding standard errors.
#'
#' @details
#' This function implements the DerSimonian-Laird method for estimating
#' heterogeneity in meta-analysis, adapted for Mendelian randomization studies.
#' The causal effect estimates (beta) and their variances are calculated using
#' the ratio method as described in Burgess and Thompson (2021), equation (4.5).
#'
#' The method-of-moments estimator is based on Cochran's Q statistic and is
#' truncated at zero to ensure non-negative variance estimates.
#'
#' @references
#' Burgess S, Thompson SG (2021). Mendelian Randomization: Methods for Causal
#' Inference Using Genetic Variants, 2nd edition. Chapman and Hall/CRC.
#'
#' @examples
#' # Example data for 5 genetic variants
#' example_data <- matrix(c(
#'   0.15, 0.10, 0.02, 0.03,
#'   0.20, 0.15, 0.03, 0.04,
#'   0.12, 0.09, 0.02, 0.03,
#'   0.18, 0.12, 0.03, 0.04,
#'   0.25, 0.18, 0.04, 0.05
#' ), ncol = 4, byrow = TRUE)
#'
#' tau2_dl(example_data)
#'
#' @export
tau2_dl <- function(data, secondorder = TRUE, beta_max = NULL, se_max = NULL) {
  gammaj_hat <- data[, 1]
  Gammaj_hat <- data[, 2]
  sigmaj_X <- data[, 3]
  sigmaj_Y <- data[, 4]

  # from Burgess and Thompson, Mendelian Randomization, 2nd edition (2021) - equation (4.5)
  betaj_hat <- Gammaj_hat/gammaj_hat
  var_betaj_hat <- (sigmaj_Y/gammaj_hat)^2
  if (secondorder) 
    var_betaj_hat <- var_betaj_hat + ((Gammaj_hat*sigmaj_X)/gammaj_hat^2)^2
  sej_hat <- sqrt(var_betaj_hat)

  beta_idx <- se_idx <- numeric(0)
  if (!is.null(beta_max))
    beta_idx <- which(abs(betaj_hat) > beta_max)
  if (!is.null(se_max))
    se_idx <- which(abs(sej_hat) > se_max)
  idx <- sort(union(beta_idx, se_idx))
  if (length(idx) > 0) {
    betaj_hat <- betaj_hat[-idx]
    var_betaj_hat <- var_betaj_hat[-idx]
    sej_hat <- sqrt(var_betaj_hat)
  }

  p <- length(betaj_hat)
  w <- 1/var_betaj_hat
  sum_w <- sum(w)
  mu_hat <- sum(betaj_hat * w)/sum_w
  W <- sum_w - (sum(w^2)/sum_w)
  Q <- sum(w*(betaj_hat - mu_hat)^2)

  tau2_hat <- max(0, (Q - (p - 1))/W)

  list(tau2_hat = tau2_hat, betaj_hat = betaj_hat, sej_hat = sej_hat)
}
