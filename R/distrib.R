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

    if(sigma_xx <= 0 || sigma_yy <= 0) {
        if(log) return(-Inf) else return(0)
    }
    
    # check positive definiteness
    det_sigma <- sigma_xx*sigma_yy - sigma_xy^2
    if (det_sigma <= 0) {
        if (log) return(-Inf) else return(0)
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
    
    if (log) return(log_dens) else return(exp(log_dens))
}
