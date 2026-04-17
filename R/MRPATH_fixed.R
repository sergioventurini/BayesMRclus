MRPATH_optimizeInitVals_fixed <- function(K, data, Nreps = 10, verbose = FALSE,
                                          altModel = FALSE, init_seed = 8686, ...) {
  initVals.list <- list()
  fit.list <- list()
  Q.vec <- rep(NA_real_, Nreps)
  
  init_m_X <- mean(data$beta.exposure)
  init_lambdaX <- sd(data$beta.exposure)
  initPis <- rep(1 / K, K)
  
  quantiles <- c(0.025, sapply(1:K, function(k) 0.025 + (0.95 / ((K - k) + 1))))
  ratio.quantiles <- quantile(
    (data$se.exposure / data$se.outcome) * (data$beta.outcome / data$beta.exposure),
    quantiles
  )
  
  for (i in seq_len(Nreps)) {
    set.seed(init_seed * i)
    initMus <- numeric(K)
    for (k in seq_len(K)) {
      initMus[k] <- runif(1, ratio.quantiles[k], ratio.quantiles[k + 1])
    }
    
    initVals <- list(
      m_X = init_m_X,
      lambdaX = init_lambdaX,
      pis = initPis,
      mus = initMus
    )
    
    if (!altModel) {
      initVals$sds <- sapply(seq_len(K), function(k)
        abs(ratio.quantiles[k + 1] - ratio.quantiles[k]) / 2)
    }
    
    fit_i <- tryCatch({
      if (altModel) {
        MRPATH::MR_PATHalt(data, initVals, ...)
      } else {
        MRPATH::MR_PATH(K, data, initVals, computeSE = FALSE, ...)
      }
    }, error = function(e) NULL)
    
    if (!is.null(fit_i)) {
      Q_i <- if (altModel) fit_i$completeDataLogLik else fit_i$convergenceInfo$completeDataLogLik
      if (is.finite(Q_i)) {
        Q.vec[i] <- Q_i
        fit.list[[i]] <- fit_i
        initVals.list[[i]] <- initVals
      }
    }
    
    if (verbose) {
      message("Run ", i, ": Q = ", Q.vec[i])
    }
  }
  
  if (all(!is.finite(Q.vec))) {
    stop("All runs failed or produced non-finite complete-data log-likelihood.")
  }
  
  best <- which.max(Q.vec)
  list(
    fit = fit.list[[best]],
    initVals = initVals.list[[best]],
    Q.vec = Q.vec,
    best = best
  )
}

#' @export
MRPATH_BIC <- function(data, K, fit) {
  p <- nrow(data)
  (-2 * fit$convergenceInfo$completeDataLogLik) + (3 * K * log(p))
}
