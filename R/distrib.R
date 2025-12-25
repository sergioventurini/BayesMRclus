dinvgamma <- function(x, alpha, beta = 1, log = FALSE) {
	if ((alpha <= 0) | (beta <= 0)) {
		stop("alpha (shape) and/or beta (scale) parameter negative in dinvgamma().\n")
	}
	log_density <- alpha*log(beta) - lgamma(alpha) - (alpha + 1)*log(x) - (beta/x)
	if (log) return(log_density) else return(exp(log_density))
}

ddirichlet <- function(x, alpha) {
	dirichlet1 <- function(x, alpha) {
		logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
		s <- sum((alpha - 1)*log(x))
		exp(sum(s) - logD)
	}
	if (!is.matrix(x)) 
		if (is.data.frame(x))
			x <- as.matrix(x)
		else x <- t(x)
	if (!is.matrix(alpha)) 
		alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), byrow = TRUE)

	if (any(alpha <= 0))
		stop("the elements of the alpha vector must be strictly positive.")
	if (any(dim(x) != dim(alpha))) 
		stop("mismatch between dimensions of x and alpha.")

	pd <- vector(length = nrow(x))
	for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, ])
	pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
	pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
	pd
}

dbivnorm <- function(x, y, mu_x, mu_y, sigma_xx, sigma_yy, sigma_xy, log = FALSE) {
  if (length(x) != length(y))
    stop("the x and y vectors must have the same length.")

  if (sigma_xx <= 0 || sigma_yy <= 0) {
    if (log) -Inf else 0
  }
  
  # check positive definiteness
  det_sigma <- sigma_xx*sigma_yy - sigma_xy^2
  if (det_sigma <= 0) {
    if (log) -Inf else 0
  }
  
  # centre observations
  dx <- x - mu_x
  dy <- y - mu_y
  
  # quadratic form: (x - mu)'*Sigma^(-1)*(x - mu)
  # note: for 2x2 matrices, the inverse is: (1/det)*[[d, -b], [-c, a]]
  inv_det <- 1/det_sigma
  quad_form <- inv_det*(sigma_yy*dx^2 - 2*sigma_xy*dx*dy + sigma_xx*dy^2)
  
  # Log density
  log_dens <- -0.5*(log(2*pi) + log(det_sigma) + quad_form)
  
  if (log) log_dens else exp(log_dens)
}

#' Half-t Distribution Density
#'
#' Computes the probability density function of the half-\eqn{t}
#' distribution.
#'
#' @param x A numeric vector of quantiles. Typically non-negative; values
#'   less than zero are allowed but correspond to the density of the underlying
#'   Student-\eqn{t} distribution reflected at zero.
#' @param alpha A positive numeric vector giving the scale parameter.
#' @param nu A positive numeric vector giving the degrees of freedom.
#' @param log Logical; if \code{TRUE}, probabilities are returned on the
#'   log scale.
#'
#' @details
#' The half-Student-\eqn{t} distribution is defined as the distribution of
#' \eqn{|T|}, where \eqn{T} follows a Student-\eqn{t} distribution with
#' \code{nu} degrees of freedom and scale parameter \code{alpha}.
#'
#' The density is
#' \deqn{
#' f(x) = \frac{2}{\mathrm{\alpha}}
#'   \frac{\Gamma((\nu + 1)/2)}{\Gamma(\nu/2)\sqrt{\pi \nu}}
#'   \left(1 + \frac{1}{\nu}\left(\frac{x}{\mathrm{\alpha}}\right)^2\right)^{-(\nu+1)/2},
#' \qquad x \ge 0.
#' }
#'
#' @return
#' A numeric vector giving the density (or log-density) evaluated at \code{x}.
#' The length of the result is the maximum of the lengths of
#' \code{x}, \code{alpha}, and \code{nu}; shorter arguments are recycled.
#'
#' @examples
#' x <- seq(0, 100, length.out = 100)
#' y <- dhalft(x, alpha = 25, nu = 5)
#' plot(x, y, type = "l", main = "Half-t density")
#'
#' # Log-density
#' dhalft(10, alpha = 25, nu = 1, log = TRUE)
#'
#' @seealso
#' \code{\link[stats]{dt}} for the Student-\eqn{t} density.
#'
#' @export
dhalft <- function(x, alpha = 25, nu = 1, log = FALSE) {
  x <- as.vector(x)
  alpha <- as.vector(alpha)
  nu <- as.vector(nu)

  if (any(alpha <= 0)) 
    stop("The alpha (scale) parameter must be positive.")

  NN <- max(length(x), length(alpha), length(nu))
  x <- rep(x, len = NN)
  alpha <- rep(alpha, len = NN)
  nu <- rep(nu, len = NN)

  log_dens <- log(2) - log(alpha) + lgamma((nu + 1)/2) - lgamma(nu/2) - 
    0.5 * log(pi * nu) - (nu + 1)/2 * log(1 + (1/nu) * (x/alpha) * 
    (x/alpha))

  if (log) log_dens else exp(log_dens)
}

pst <- function(q, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE)  {
  q <- as.vector(q)
  mu <- as.vector(mu)
  sigma <- as.vector(sigma)
  nu <- as.vector(nu)
  if (any(sigma <= 0)) 
    stop("The sigma parameter must be positive.")
  if (any(nu <= 0)) 
    stop("The nu parameter must be positive.")

  NN <- max(length(q), length(mu), length(sigma), length(nu))
  q <- rep(q, len = NN)
  mu <- rep(mu, len = NN)
  sigma <- rep(sigma, len = NN)
  nu <- rep(nu, len = NN)
  p <- pt({
      q - mu
  }/sigma, df = nu, lower.tail = lower.tail, log.p = log.p)
  temp <- which(nu > 1e+06)
  p[temp] <- pnorm(q[temp], mu[temp], sigma[temp], lower.tail = lower.tail, 
      log.p = log.p)

  p
}

qst <- function(p, mu = 0, sigma = 1, nu = 10, lower.tail = TRUE, log.p = FALSE) {
  p <- as.vector(p)
  mu <- as.vector(mu)
  sigma <- as.vector(sigma)
  nu <- as.vector(nu)
  if (any(p < 0) || any(p > 1)) 
    stop("p must be in [0,1].")
  if (any(sigma <= 0)) 
    stop("The sigma parameter must be positive.")
  if (any(nu <= 0)) 
    stop("The nu parameter must be positive.")

  NN <- max(length(p), length(mu), length(sigma), length(nu))
  p <- rep(p, len = NN)
  mu <- rep(mu, len = NN)
  sigma <- rep(sigma, len = NN)
  nu <- rep(nu, len = NN)
  q <- mu + sigma * qt(p, df = nu, lower.tail = lower.tail)
  temp <- which(nu > 1e+06)
  q[temp] <- qnorm(p[temp], mu[temp], sigma[temp], lower.tail = lower.tail, log.p = log.p)

  q
}

qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  if (any(p < 0) || any(p > 1)) 
    stop("p must be in [0,1].")
  if (a >= b) 
    stop("Lower bound a is not less than upper bound b.")

  q <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  q <- Gin(G(a, ...) + p * {G(b, ...) - G(a, ...)}, ...)

  q
}

rtrunc <-function(n, spec, a = -Inf, b = Inf, ...) {
  if (a >= b) 
    stop("Lower bound a is not less than upper bound b.")

  x <- u <- runif(n)
  x <- qtrunc(u, spec, a = a, b = b, ...)

  x
}

#' Random Generation from the Half-t Distribution
#'
#' Generate random deviates from the half-t distribution, which is a 
#' truncated Student's t-distribution restricted to positive values.
#'
#' @param n integer; number of observations to generate.
#' @param alpha numeric; scale parameter. Must be positive. Default is 25.
#'   Vectors are recycled to length \code{n}.
#' @param nu numeric; degrees of freedom parameter. Must be positive. 
#'   Default is 1. Vectors are recycled to length \code{n}.
#'
#' @return A numeric vector of length \code{n} containing random deviates 
#'   from the half-t distribution.
#'
#' @details
#' The half-t distribution is obtained by truncating a Student's t-distribution 
#' with location parameter 0, scale parameter \code{alpha}, and \code{nu} 
#' degrees of freedom to the positive half of the real line.
#'
#' @seealso \code{\link{EXhalft}} for the expected value of the half-t distribution.
#'
#' @examples
#' # Generate 10 random values with default parameters
#' rhalft(10)
#' 
#' # Generate 100 values with custom parameters
#' x <- rhalft(100, alpha = 10, nu = 5)
#' hist(x, main = "Half-t Distribution")
#'
#' @export
rhalft <- function(n, alpha = 25, nu = 1) {
  alpha <- rep(alpha, len = n)
  nu <- rep(nu, len = n)
  if (any(alpha <= 0)) 
    stop("The alpha parameter must be positive.")
  x <- rtrunc(n, "st", a = 0, b = Inf, mu = 0, sigma = alpha, 
    nu = nu)
  
  x
}

#' Expected Value of the Half-t Distribution
#'
#' Calculate the expected value (mean) of the half-t distribution.
#'
#' @param alpha numeric; scale parameter. Must be positive. Default is 25.
#' @param nu numeric; degrees of freedom parameter. Must be greater than 1 
#'   for the expected value to be finite. Default is 1.
#'
#' @return The expected value of the half-t distribution. Returns \code{Inf} 
#'   when \code{nu <= 1}, as the expected value does not exist in this case.
#'
#' @details
#' The expected value of the half-t distribution is given by:
#' \deqn{E[X] = 2\alpha\sqrt{\frac{\nu}{\pi}} \frac{\Gamma((\nu+1)/2)}{\Gamma(\nu/2)(\nu-1)}}
#' for \eqn{\nu > 1}. When \eqn{\nu \leq 1}, the expected value is infinite.
#'
#' @seealso \code{\link{rhalft}} for random generation from the half-t distribution.
#'
#' @examples
#' # Expected value with default parameters (returns Inf since nu = 1)
#' EXhalft()
#' 
#' # Expected value with nu > 1
#' EXhalft(alpha = 10, nu = 5)
#' 
#' # Verify with simulation
#' set.seed(123)
#' mean(rhalft(10000, alpha = 10, nu = 5))
#' EXhalft(alpha = 10, nu = 5)
#'
#' @export
EXhalft <- function(alpha = 25, nu = 1) {
  if (nu <= 1) {
    Inf
  }
  else {
    2*alpha*sqrt(nu/pi)*gamma((nu + 1)/2)/(gamma(nu/2)*(nu - 1))
  }
}
