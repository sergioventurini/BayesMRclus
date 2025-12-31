#' Function to compute the starting values before fitting a BayesMR models.
#'
#' \code{bayesmr_init()} is the main function that estimates a BayesMR model.
#'
#' @param data A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param p A length-one numeric vector indicating the number of dimensions of the
#'   latent space.
#' @param G A length-one numeric vector indicating the number of cluster to
#'   partition the \emph{S} subjects.
#' @param random.start A length-one logical vector. If \code{TRUE} the starting
#'   values are drawn randomly, otherwise.
#' @param partition A length-one numeric vector providing the user-defined
#'   starting partition.
#' @return A named \code{list} with the following items:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates starting values}
#'     \item{\code{x}: }{numeric vector of initial cluster memberships}
#'     \item{\code{ng}: }{numeric vector of initial cluster sizes}
#'     \item{\code{alpha}: }{numeric vector of alpha starting values}
#'     \item{\code{eta}: }{numeric vector of eta starting values}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 starting values}
#'     \item{\code{lambda}: }{numeric vector of lambda starting values}
#'   }
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bayesmr}()} for fitting a BayesMR model.
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#' @examples
#' data(simdiss, package = "bayesmr")
#' bayesmr_init(simdiss@diss, p = 2, G = 3, random.start = TRUE)
#' @export
bayesmr_init <- function(data, p, G, random.start, partition) {
  n <- nrow(data)

  # # initialize x (cluster labels)
  # if (random.start) {
  #   x <- sample(1:G, S, replace = TRUE)
  #   while (length(unique(x)) < G) {
  #     x <- sample(1:G, S, replace = TRUE)
  #   }
  # } else {
  #   if (is.null(partition)) {
  #     dmat <- list2matrix(data)
  #     d.clust <- hclust(dist(dmat, method = method), method = "ward.D")
  #     x <- cutree(d.clust, k = G)
  #   } else {
  #     if (length(partition) != S) {
  #       stop(paste0("the initial partition must include S = ", S, " values."))
  #     }
  #     if (length(unique(partition)) != G) {
  #       stop(paste0("the initial partition must include G = ", G, " unique values."))
  #     }
  #     x <- partition
  #   }
  # }
  # ng <- table(factor(x, levels = 1:G))
  
  # Dm <- list2array(D)
  # z <- array(NA, dim = c(n, p, G))
  # alpha <- eta <- sigma2 <- dsigma2 <- numeric(G)

  # for (g in 1:G) {
  #   # initialize Z_g
  #   Dg <- Dm[, , x == g]
  #   Dm_bar <- apply(Dg, c(1, 2), mean)
  #   D_bar <- as.dist(Dm_bar)
  #   d_ov.mean <- mean(D_bar)
  #   Dm_above <- as.matrix(as.dist(Dm_bar > d_ov.mean))
    
  #   Dg <- Dm[, , x == g]
  #   Dm_sum <- apply(Dg, c(1, 2), sum)
  #   z_mds <- stats::cmdscale(d = Dm_sum, k = p)
  #   if (ncol(z_mds) != p) {
  #     stop(paste0("the initialization of the MDS configuration for cluster g = ", g, " of ", G,
  #       " using p = ", p, " dimensions failed."))
  #   }
  #   z[, , g] <- scale(z_mds)
    
  #   # initialize alpha_g
  #   alpha.glm <- glm(as.numeric(as.dist(Dm_above)) ~ 1, family = family, offset = as.numeric(dist(z[, , g])))
  #   alpha[g] <- coef(alpha.glm)
    
  #   # initialize sigma2_g
  #   sigma2[g] <- vcov(alpha.glm)[1, 1]
  # }
  
  # initialize psi and tau
  psi_tau <- rhalft(2, alpha = .01, nu = 3)
  psi <- psi_tau[1]
  tau <- psi_tau[2]

  # initialize gamma
  gamma <- rnorm(1, mean = 0, sd = .1)

  # initialize beta
  beta <- rnorm(1, mean = 0, sd = .1)
  
  list(gamma = gamma, beta = beta, psi = psi, tau = tau)
}
