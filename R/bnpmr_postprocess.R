# =============================================================================
# 1. Read MCMC output
# =============================================================================

#' Read all CSV output files written by mr_mixture.cpp
#'
#' @param res        Output of the bayesmr_mix_het() function
#' @param burnin     Number of initial iterations to discard
#' @param pleiotropy Declare whether the model assumes pleiotropy
#' @return           A list with components described in the header above
read_mcmc_output <- function(res, burnin, pleiotropy = TRUE) {
  nchains <- length(res@results)
  
  gamma <- psi <- alpha <- numeric(0)
  if (pleiotropy) tau <- numeric(0)
  chain_list <- list()
  for (c in 1:nchains) {
    gamma <- res@results[[c]]@gamma.chain
    psi   <- res@results[[c]]@psi.chain
    if (pleiotropy) tau <- res@results[[c]]@tau.chain
    alpha <- res@results[[c]]@alpha.chain
  
    S <- length(gamma)
    if (pleiotropy) {
      stopifnot(length(psi) == S, length(tau) == S, length(alpha) == S)
    } else {
      stopifnot(length(psi) == S, length(alpha) == S)
    }
    
    # xi: S × p integer matrix
    xi_raw <- res@results[[c]]@xi.chain
    xi <- as.matrix(xi_raw)
    storage.mode(xi) <- "integer"
    p <- ncol(xi)
    stopifnot(nrow(xi) == S)
    
    # beta: variable-length rows, NA-padded
    beta_raw <- res@results[[c]]@beta.chain

    K_max <- max(apply(beta_raw, 1, function(x) sum(!is.na(x))), na.rm = TRUE)
    beta_raw <- beta_raw[, 1:K_max]
    beta_mat <- as.matrix(beta_raw)
    storage.mode(beta_mat) <- "double"
    stopifnot(nrow(beta_mat) == S)
    
    # Convert each row to a trimmed vector (drop trailing NAs)
    beta <- lapply(seq_len(S), function(s) {
      b <- beta_mat[s, ]
      b[!is.na(b)]
    })
    
    if (pleiotropy) {
      chain_list[[c]] <- list(gamma = gamma, psi = psi, tau = tau, alpha = alpha,
           xi = xi, beta = beta, S = S, p = p)
    } else {
      chain_list[[c]] <- list(gamma = gamma, psi = psi, alpha = alpha,
           xi = xi, beta = beta, S = S, p = p)
    }

    # remove burnin
    chain_list[[c]] <- remove_burnin(chain_list[[c]], burnin, pleiotropy = pleiotropy)
  }
  combined <- combine_chains(chain_list, pleiotropy = pleiotropy)
  
  list(combined = combined, chains = chain_list)
}

#' Combine a list of chain objects into one
#'
#' All chains must have the same p.  xi labels are chain-local and remain
#' so after concatenation — this is harmless because all post-processing
#' functions use only within-iteration label consistency.
#'
#' @param chain_list  List of out objects from remove_burnin()
#' @param pleiotropy  Declare whether the model assumes pleiotropy
#' @return            A single merged out object
combine_chains <- function(chain_list, pleiotropy = TRUE) {
  p <- chain_list[[1]]$p
  
  out <- list(
    gamma = do.call(c,     lapply(chain_list, `[[`, "gamma")),
    psi   = do.call(c,     lapply(chain_list, `[[`, "psi")),
    alpha = do.call(c,     lapply(chain_list, `[[`, "alpha")),
    xi    = do.call(rbind, lapply(chain_list, `[[`, "xi")),
    beta  = do.call(c,     lapply(chain_list, `[[`, "beta")),
    S     = sum(sapply(chain_list, `[[`, "S")),
    p     = p
  )

  if (pleiotropy)
    out <- append(out, list(tau = do.call(c, lapply(chain_list, `[[`, "tau"))), after = 2)

  out
}


# =============================================================================
# 2. Burn-in removal
# =============================================================================

#' Discard the first `burnin` iterations
#'
#' @param out        Object from which to remove the burnin
#' @param burnin     Number of initial iterations to discard
#' @param pleiotropy Declare whether the model assumes pleiotropy
#' @return           Trimmed object with updated S
remove_burnin <- function(out, burnin, pleiotropy = TRUE) {
  stopifnot(burnin < out$S)
  idx <- seq(burnin + 1L, out$S)
  out$gamma <- out$gamma[idx]
  out$psi   <- out$psi[idx]
  if (pleiotropy) out$tau <- out$tau[idx]
  out$alpha <- out$alpha[idx]
  out$xi    <- out$xi[idx, , drop = FALSE]
  out$beta  <- out$beta[idx]
  out$S     <- length(idx)
  out
}


# =============================================================================
# 3. Convergence diagnostics
# =============================================================================

#' Basic convergence diagnostics for scalar parameters
#'
#' Returns a data frame with posterior mean, SD, and the Gelman-Rubin Rhat
#' if multiple chains are provided (single chain: NA).
#' Also returns ESS via the batch-means estimator.
#'
#' @param out        Object from remove_burnin() (single chain)
#' @param chains     Optional list of additional out objects (for multi-chain Rhat)
#' @param pleiotropy Declare whether the model assumes pleiotropy
#' @return           Data frame with columns: param, mean, sd, ess, rhat
convergence_diagnostics <- function(out, chains = NULL, pleiotropy = TRUE) {
  
  # Effective sample size via batch means
  ess_bm <- function(x, b = max(1L, floor(length(x)^(1/3)))) {
    n  <- length(x); nb <- floor(n / b)
    bm <- sapply(seq_len(nb), function(i) mean(x[((i-1)*b+1):(i*b)]))
    n * var(x) / (b * var(bm) + 1e-15)
  }
  
  K_vec <- sapply(out$beta, length)
  
  params <- list(
    gamma = out$gamma,
    psi   = out$psi,
    alpha = out$alpha,
    K     = K_vec
  )

  if (pleiotropy)
    params <- append(params, list(tau = out$tau), after = 2)
  
  rows <- lapply(names(params), function(nm) {
    x <- params[[nm]]
    data.frame(param = nm,
               mean  = mean(x),
               sd    = sd(x),
               ess   = round(ess_bm(x)),
               rhat  = NA_real_,
               stringsAsFactors = FALSE)
  })
  diag_df <- do.call(rbind, rows)
  
  # Multi-chain Rhat (Gelman-Rubin)
  if (!is.null(chains)) {
    all_chains <- c(list(out), chains)
    for (nm in c("gamma", "psi", "tau", "alpha")) {
      chains_list <- lapply(all_chains, function(o) get(nm, envir = as.environment(o)))
      m  <- length(chains_list)
      n  <- min(sapply(chains_list, length))
      chains_list <- lapply(chains_list, function(x) x[seq_len(n)])
      W  <- mean(sapply(chains_list, var))
      gm <- mean(sapply(chains_list, mean))
      B  <- n * var(sapply(chains_list, mean))
      Vhat <- (n - 1) / n * W + (m + 1) / (m * n) * B
      rhat <- sqrt(Vhat / (W + 1e-15))
      diag_df$rhat[diag_df$param == nm] <- round(rhat, 4)
    }
  }
  
  rownames(diag_df) <- NULL
  diag_df
}


# =============================================================================
# 4. Posterior distribution of K
# =============================================================================

#' Posterior distribution of the number of active clusters K
#'
#' @param out   Object from remove_burnin()
#' @return      Data frame with columns k, count, proportion
posterior_K <- function(out) {
  K_vec <- sapply(out$beta, length)
  tab   <- table(K_vec)
  df    <- data.frame(k           = as.integer(names(tab)),
                      count       = as.integer(tab),
                      proportion  = as.numeric(tab) / out$S,
                      stringsAsFactors = FALSE)
  df[order(df$k), ]
}


# =============================================================================
# 5. Posterior Similarity Matrix
# =============================================================================

#' Compute the p × p Posterior Similarity Matrix
#'
#' \eqn{PSM_{i,j}} = proportion of iterations in which SNPs i and j
#'                  are assigned to the same cluster.
#'
#' @param out   Object from remove_burnin()
#' @return      Symmetric numeric matrix, dim p × p, entries in [0,1]
compute_psm <- function(out) {
  S <- out$S; p <- out$p
  PSM <- matrix(0.0, p, p)
  for (s in seq_len(S)) {
    xi_s <- out$xi[s, ]
    # For each cluster label, add 1 to all pairs in that cluster
    for (k in unique(xi_s)) {
      idx <- which(xi_s == k)
      PSM[idx, idx] <- PSM[idx, idx] + 1.0
    }
  }
  PSM / S
}


# =============================================================================
# 6. Representative partition via Binder loss minimisation
# =============================================================================

#' Find the MCMC sample that minimises Binder's expected loss
#' relative to the PSM
#'
#' Binder loss for a partition c is proportional to
#'   \eqn{\sum_{i<j} { PSM_{i,j} \cdot 1(c_i \ne c_j)  +
#'               (1 - PSM_{i,j}) \cdot 1(c_i = c_j) }}
#'
#' @param out   Object from remove_burnin()
#' @param loss  salso loss function; default is VI() as in Dahl et al. (2022)
#'              Use binder() for Binder loss with automatic a optimization
#' @return      List with elements:
#'                $partition  integer vector length p
#'                $loss       expected loss of the selected partition
#'                $info       salso info data frame
salso_partition <- function(out, loss = VI()) {
  # salso expects an integer matrix with rows = iterations, cols = items
  # out$xi is already in this format (S × p, 1-indexed)
  est  <- salso::salso(out$xi, loss = loss)
  info <- attr(est, "info")
  list(
    partition = as.integer(est),
    loss      = info$expectedLoss,
    info      = info
  )
}


# =============================================================================
# 7. Per-SNP BMA causal effects
# =============================================================================

#' Bayesian model averaging posterior for each SNP's causal effect
#'
#' At each iteration s, beta_j^(s) = beta*_{xi_j^(s)}.
#' Averages these over all post-burn-in samples.
#'
#' @param out   Object from remove_burnin()
#' @param ci    Credible interval level (default 0.95)
#' @return      Data frame with columns:
#'                snp_index, mean, sd, lower, upper,
#'                prob_positive, prob_negative
snp_bma_effects <- function(out, ci = 0.95) {
  S <- out$S; p <- out$p
  alpha_ci <- (1 - ci) / 2
  
  # Matrix S × p of SNP-specific causal effects
  beta_snp <- matrix(NA_real_, S, p)
  for (s in seq_len(S)) {
    beta_snp[s, ] <- out$beta[[s]][out$xi[s, ]]
  }
  
  data.frame(
    snp_index        = seq_len(p),
    mean             = colMeans(beta_snp),
    sd               = apply(beta_snp, 2, sd),
    lower            = apply(beta_snp, 2, quantile, probs = alpha_ci),
    upper            = apply(beta_snp, 2, quantile, probs = 1 - alpha_ci),
    prob_positive    = colMeans(beta_snp > 0),
    prob_negative    = colMeans(beta_snp < 0),
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# 8. Cluster-specific effect summaries
# =============================================================================

#' Summarise cluster-specific effects from the Binder partition
#'
#' For each cluster in the representative partition, reports the
#' posterior mean, SD, and credible interval of the cluster atom beta*_k,
#' marginalised over the label-switching uncertainty by using the
#' per-SNP BMA effects restricted to SNPs assigned to that cluster.
#'
#' @param binder  Output from salso_partition()
#' @param bma     Output from snp_bma_effects()
#' @param out     Output from remove_burnin()
#' @param ci      Credible interval level (default 0.95)
#' @return        Data frame with one row per cluster
cluster_summaries <- function(binder, bma, out, ci = 0.95) {
  part   <- binder$partition
  K_rep  <- length(unique(part))
  labs   <- sort(unique(part))
  alpha_ci <- (1 - ci) / 2
  
  # Per-SNP draws for the Binder iteration
  beta_snp_all <- matrix(NA_real_, out$S, out$p)
  for (s in seq_len(out$S))
    beta_snp_all[s, ] <- out$beta[[s]][out$xi[s, ]]
  
  rows <- lapply(seq_along(labs), function(i) {
    k      <- labs[i]
    snps_k <- which(part == k)
    # Pool all per-SNP draws within this cluster
    pool   <- as.vector(beta_snp_all[, snps_k])
    data.frame(
      cluster   = i,
      n_snps    = length(snps_k),
      mean      = mean(pool),
      sd        = sd(pool),
      lower     = quantile(pool, alpha_ci),
      upper     = quantile(pool, 1 - alpha_ci),
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)
  rownames(rows) <- NULL
  rows
}


# =============================================================================
# 9. Pleiotropy check: dispersion statistic T
# =============================================================================

#' Posterior predictive dispersion statistic for pleiotropy assessment
#'
#' For each posterior draw (beta_1,...,beta_p, tau^2), computes the
#' standardised residuals r_j = (Gamma_hat_j - beta_j * gamma_hat_j) /
#' sqrt(s^2_Gamma_j + tau^2) and the dispersion statistic
#' T^(s) = (1/p) * sum_j r_j^2.
#' T^(s) should concentrate near 1 under correct model specification.
#'
#' @param out     Object from remove_burnin()
#' @param dat     Data frame with columns Gamma_hat, gamma_hat,
#'                sigma2Y (= s^2_Gamma_j, i.e. variance not SE)
#' @return        List with:
#'                  $T_draws  numeric vector length S of T^(s) values
#'                  $mean_T   posterior mean of T
#'                  $ci_T     95% credible interval for T
#'                  $r_matrix S × p matrix of standardised residuals
pleiotropy_check <- function(out, dat) {
  S  <- out$S
  p  <- out$p
  Gh <- dat$Gamma_hat
  gh <- dat$gamma_hat
  s2Y<- dat$sigma2Y
  
  # Per-SNP BMA beta draws (S × p)
  beta_snp <- matrix(NA_real_, S, p)
  for (s in seq_len(S))
    beta_snp[s, ] <- out$beta[[s]][out$xi[s, ]]
  
  T_vec    <- numeric(S)
  r_matrix <- matrix(NA_real_, S, p)
  
  for (s in seq_len(S)) {
    if (!is.null(out$tau)) {
      tau2 <- out$tau[s]^2
    } else {
      tau2 <- 0
    }
    fitted <- beta_snp[s, ] * gh          # beta_j * gamma_hat_j
    resid  <- Gh - fitted
    denom  <- sqrt(s2Y + tau2)
    r      <- resid / denom
    r_matrix[s, ] <- r
    T_vec[s]      <- mean(r^2)
  }
  
  list(T_draws  = T_vec,
       mean_T   = mean(T_vec),
       ci_T     = quantile(T_vec, c(0.025, 0.975)),
       r_matrix = r_matrix)
}


# =============================================================================
# 10. Figures
# =============================================================================

#' Plot the Posterior Similarity Matrix as a heatmap
#'
#' Rows and columns are reordered by the Binder partition (SNPs in the same
#' cluster are grouped together), and a colour strip is added on both margins
#' indicating cluster membership.
#'
#' @param psm        Matrix from compute_psm()
#' @param binder     Output from salso_partition()
#' @param snp_names  Optional character vector of SNP names (length p)
#' @param name_cex   Character expansion for SNP names
#' @param file       Optional file path; if NULL the plot goes to the current device
#' @param width      Width of the image file saved on disk
#' @param height     Height of the image file saved on disk
#' @param res        Resolution of the image file saved on disk
plot_psm <- function(psm, binder, snp_names = NULL,
                     name_cex     = 0.35,
                     file         = NULL,
                     width        = 1600,
                     height       = 1500,
                     res          = 220) {

  # ── Input checks ────────────────────────────────────────────────────────────
  p <- nrow(psm)
  stopifnot(ncol(psm) == p, length(binder$partition) == p)
  if (!is.null(snp_names)) stopifnot(length(snp_names) == p)

  # ── Reorder by partition ─────────────────────────────────────────────────
  part    <- binder$partition
  ord     <- order(part)
  psm_ord <- psm[ord, ord]
  K       <- length(unique(part))

  # SNP labels in reordered order (fall back to index if no names supplied)
  labs <- if (!is.null(snp_names)) snp_names[ord] else as.character(ord)

  # ── Colour palette ───────────────────────────────────────────────────────
  strip_cols <- c("#D62728","#1F77B4","#2CA02C","#FF7F0E",
                  "#9467BD","#8C564B","#E377C2","#7F7F7F")
  strip_cols <- rep(strip_cols, length.out = K)
  snp_col    <- strip_cols[part[ord]]

  # ── Margins: wider when names are shown ─────────────────────────────────
  # The colour strip adds ~4% of p units on each margin side; we also need
  # room for rotated text.  We use large bottom and left margins and let
  # axis() handle the tick labels.
  mar_base <- if (!is.null(snp_names)) c(5, 5, 2, 2) else c(5, 5, 2, 2)

  if (!is.null(file)) png(file, width = width, height = height, res = res)
  op <- par(mar = mar_base)

  # ── Main heatmap ─────────────────────────────────────────────────────────
  image(seq_len(p), seq_len(p),
        psm_ord,
        col  = colorRampPalette(c("white", "#08306b"))(100),
        axes = FALSE,
        xlab = "",
        ylab = "",
        # main = "Posterior Similarity Matrix")
        main = "")

  # ── Axes ──────────────────────────────────────────────────────────────────
  if (!is.null(snp_names)) {
    # x-axis: SNP names rotated 90 degrees (srt=90, las=2)
    axis(1, at = seq_len(p), labels = labs,
         las = 2, cex.axis = name_cex, tick = FALSE, gap.axis = 0)

    # y-axis: SNP names horizontal (las=2)
    axis(2, at = seq_len(p), labels = labs,
         las = 2, cex.axis = name_cex, tick = FALSE, gap.axis = 0)
  } else {
    # Fall back to numeric indices
    axis(1, at = pretty(seq_len(p)), cex.axis = 0.7)
    axis(2, at = pretty(seq_len(p)), cex.axis = 0.7)
  }
  title(xlab = "SNP index (reordered)", ylab = "SNP index (reordered)")

  # ── Colour strip on margins (outside the plot region) ────────────────────
  # Strip width is expressed in data units (fraction of p)
  strip_w <- p * 0.03          # width of the colour strip in data units
  strip_gap <- p * 0.005       # small gap between strip and plot edge

  for (i in seq_len(p)) {
    # Bottom strip (below x-axis)
    rect(i - 0.5, -(strip_gap + strip_w), i + 0.5, -strip_gap,
         col = snp_col[i], border = NA, xpd = TRUE)
    # Left strip (left of y-axis)
    rect(-(strip_gap + strip_w), i - 0.5, -strip_gap, i + 0.5,
         col = snp_col[i], border = NA, xpd = TRUE)
  }

  # ── Cluster boundary lines ───────────────────────────────────────────────
  counts  <- rle(part[ord])$lengths          # rle on already-sorted part
  cum_cnt <- cumsum(counts)
  for (b in head(cum_cnt, -1))
    abline(h = b + 0.5, v = b + 0.5, col = "white", lwd = 1.2)

  # ── Legend ───────────────────────────────────────────────────────────────
  # bty = "n" (not "0" which is not a valid value and falls back silently)
  legend("topleft",
         legend = paste("Cluster", seq_len(K)),
         fill   = strip_cols[seq_len(K)],
         bty    = "o",
         bg     = "white",
         cex    = 0.85)

  par(op)
  if (!is.null(file)) { dev.off(); invisible(file) }
}


#' Plot per-SNP BMA posterior means with credible intervals
#'
#' SNPs are ordered by their Wald ratio and coloured by their Binder
#' partition cluster.  The empirical Wald ratio is overlaid.
#'
#' @param bma       Output from snp_bma_effects()
#' @param binder    Output from salso_partition()
#' @param snp_names Optional character vector of SNP names (length p)
#' @param wald      Numeric vector of Wald ratios (length p)
#' @param ordering  Ordering of the causal effects (default "post_mean")
#' @param name_cex   Character expansion for SNP names
#' @param file       Optional file path; if NULL the plot goes to the current device
#' @param width      Width of the image file saved on disk
#' @param height     Height of the image file saved on disk
#' @param res        Resolution of the image file saved on disk
plot_snp_effects <- function(bma, binder,
                             snp_names  = NULL,
                             wald,
                             ordering   = "post_mean",
                             name_cex   = 0.35,
                             file       = NULL,
                             width      = 1600,
                             height     = 1100,
                             res        = 220) {

  # ── Input checks ────────────────────────────────────────────────────────────
  p <- nrow(bma)
  stopifnot(length(binder$partition) == p, length(wald) == p)
  if (!is.null(snp_names)) stopifnot(length(snp_names) == p)

  # ── Ordering ─────────────────────────────────────────────────────────────
  ord <- switch(ordering,
    wald      = order(wald),
    post_mean = order(bma$mean),
    seq_len(p)          # default: original order
  )
  if (is.null(ord))
    stop("'ordering' must be \"wald\", \"post_mean\", or \"none\".")

  # SNP labels at each rank position (reordered)
  # Position i on the x-axis corresponds to snp_names[ord[i]]
  x_labs <- if (!is.null(snp_names)) snp_names[ord] else as.character(ord)

  # ── Cluster colours ───────────────────────────────────────────────────────
  part       <- binder$partition
  K          <- length(unique(part))
  strip_cols <- c("#D62728","#1F77B4","#2CA02C","#FF7F0E",
                  "#9467BD","#8C564B","#E377C2","#7F7F7F")
  strip_cols <- rep(strip_cols, length.out = K)
  pt_col     <- strip_cols[part]   # colour for each SNP by its cluster

  # ── Y-axis range ──────────────────────────────────────────────────────────
  ylim <- range(c(bma$lower, bma$upper, wald), na.rm = TRUE) * 1.05

  # ── Margins: extra bottom room when names are shown ──────────────────────
  # With las=2, x-axis labels are rotated 90° and extend downward.
  # Increase bottom margin proportionally to the longest name.
  if (!is.null(snp_names)) {
    max_nchar <- max(nchar(x_labs))
    # Each character is roughly 0.07 lines at cex=1; scale by name_cex
    bot_mar <- max(5, ceiling(max_nchar * 0.07 / name_cex))
  } else {
    bot_mar <- 5
  }

  if (!is.null(file)) png(file, width = width, height = height, res = res)
  op <- par(mar = c(bot_mar, 5, 2, 2))

  # ── Base plot (no x-axis: we draw it manually below) ─────────────────────
  xlab_str <- if (!is.null(snp_names)) {
    ""    # axis label replaced by individual SNP names
  } else {
    paste0("SNP rank (by ",
           ifelse(ordering == "wald", "Wald ratio", "posterior mean"), ")")
  }

  plot(NA,
       xlim = c(0.5, p + 0.5), ylim = ylim,
       xlab = xlab_str,
       ylab = "Causal effect",
       # main = "Per-SNP BMA posterior causal effects",
       main = "",
       xaxt = "n")   # suppress default x-axis
  abline(h = 0, col = "grey60", lty = 2)

  # ── X-axis ────────────────────────────────────────────────────────────────
  if (!is.null(snp_names)) {
    # Individual SNP names, rotated 90 degrees, at each rank position
    op2 <- par(mgp = c(3, 0.3, 0))
    axis(1,
         at       = seq_len(p),
         labels   = x_labs,
         las      = 2,           # perpendicular to axis => vertical text
         cex.axis = name_cex,
         tick     = FALSE,
         gap.axis = 0)
    par(op2)
  } else {
    # Plain numeric axis (rank numbers)
    axis(1)
  }
  title(xlab = paste0("SNP rank (by ",
                      ifelse(ordering == "wald",
                             "Wald ratio", "posterior mean"), ")"))

  # ── Data layers ───────────────────────────────────────────────────────────
  # Background: Wald ratios
  points(seq_len(p), wald[ord], pch = 1, cex = 0.3, col = "grey70")

  # Credible interval segments
  for (i in seq_len(p)) {
    j <- ord[i]
    segments(i, bma$lower[j], i, bma$upper[j],
             col = pt_col[j], lwd = 0.8)
  }

  # Posterior mean points
  points(seq_len(p), bma$mean[ord], pch = 19, cex = 0.4, col = pt_col[ord])

  # ── Legend ────────────────────────────────────────────────────────────────
  legend("topleft",
         legend = c(paste("Cluster", seq_len(K)), "Wald ratio"),
         col    = c(strip_cols[seq_len(K)], "grey70"),
         pch    = c(rep(19, K), 1),
         bty    = "n",
         cex    = 0.65)

  par(op)
  if (!is.null(file)) { dev.off(); invisible(file) }
}


#' Plot posterior distribution of K
#'
#' @param Kd    Output from posterior_K()
#' @param file  Optional file path
plot_posterior_K <- function(Kd, file = NULL) {
  if (!is.null(file)) png(file, width = 1200, height = 900, res = 220)
  
  par(mar = c(5, 5, 2, 2))
  barplot(Kd$proportion,
          names.arg = Kd$k,
          col  = "#1F77B4",
          xlab = "Number of clusters K",
          ylab = "Posterior probability",
          # main = "Posterior distribution of K",
          main = "",
          ylim = c(0, max(Kd$proportion) * 1.15))
  box()
  
  if (!is.null(file)) { dev.off(); invisible(file) }
}


#' Trace plots for scalar parameters
#'
#' @param out        Object from remove_burnin()
#' @param pleiotropy Declare whether the model assumes pleiotropy
#' @param file       Optional file path
plot_traces <- function(out, file = NULL, pleiotropy = TRUE) {
  if (!is.null(file)) png(file, width = 1800, height = 1400, res = 220)
  
  K_vec <- sapply(out$beta, length)
  if (pleiotropy) {
    op <- par(mfrow = c(3, 2), mar = c(4.5, 4, 3, 1))
  } else {
    op <- par(mfrow = c(2, 2), mar = c(4.5, 4, 3, 1))
  }
  
  for (nm in list(
    list(x = out$gamma, lab = expression(gamma)),
    list(x = out$psi,   lab = expression(psi)),
    list(x = out$alpha, lab = expression(alpha)),
    list(x = K_vec,     lab = expression("K"))
  )) {
    ttl <- if (is.expression(nm$lab)) {
      bquote("Trace:" ~ .(nm$lab[[1]]))
    } else {
      paste("Trace:", nm$lab)
    }
    plot(nm$x, type = "l", col = "#1F77B4", lwd = 0.4,
         xlab = "Iteration", ylab = nm$lab, main = ttl)
    abline(h = mean(nm$x), col = "red", lty = 2, lwd = 1.2)
  }
  
  if (pleiotropy) {
    nm <- list(x = out$tau,   lab = expression(tau))
    ttl <- if (is.expression(nm$lab)) {
      bquote("Trace:" ~ .(nm$lab[[1]]))
    } else {
      paste("Trace:", nm$lab)
    }
    plot(nm$x, type = "l", col = "#1F77B4", lwd = 0.4,
         xlab = "Iteration", ylab = nm$lab, main = ttl)
    abline(h = mean(nm$x), col = "red", lty = 2, lwd = 1.2)

    # ACF for tau
    acf(out$tau, main = expression("ACF: " * tau),
        col = "#1F77B4", lwd = 1.5)
  }
  
  par(op)
  if (!is.null(file)) { dev.off(); invisible(file) }
}


# =============================================================================
# 11. Full pipeline wrapper
# =============================================================================

#' Run the complete post-processing pipeline
#'
#' Reads output, removes burn-in, computes all summaries, and
#' optionally produces figures.
#'
#' @param res        Output of the bayesmr_mix_het() function
#' @param burnin     Number of burn-in iterations to discard
#' @param dat        Data frame with columns Gamma_hat, gamma_hat, sigma2Y, sigma2X
#'                   (needed for pleiotropy_check; pass NULL to skip)
#' @param fig_dir    Directory for output figures (NULL = no figures)
#' @param ci         Credible interval level (default 0.95)
#' @param ordering   Ordering of the causal effects (default "post_mean")
#' @param loss       Loss to use to identify optimal partition (default "binder")
#' @param verbose    Display the process (default TRUE)
#' @param snp_names  Optional character vector of SNP names (length p)
#' @param name_cex   Character expansion for SNP names
#' @param pleiotropy Declare whether the model assumes pleiotropy
#' @return           Named list with all inferential summaries
#' @export
run_bnpmr_postprocess <- function(res,
                                  burnin     = 2000L,
                                  dat        = NULL,
                                  fig_dir    = NULL,
                                  ci         = 0.95,
                                  ordering   = "post_mean",
                                  loss       = "binder",
                                  verbose    = TRUE,
                                  snp_names  = NULL,
                                  name_cex   = 0.35,
                                  pleiotropy = TRUE) {
  if (verbose) cat("Reading MCMC output ...\n")
  out  <- read_mcmc_output(res, burnin, pleiotropy = pleiotropy)
  if (verbose) cat(sprintf("  S = %d iterations,  p = %d SNPs\n", out$S, out$p))
  
  if (verbose) cat("Convergence diagnostics ...\n")
  diag <- lapply(out$chain_list, convergence_diagnostics, pleiotropy = pleiotropy)
  if (verbose) print(diag)
  
  out <- out$combined
  
  if (verbose) cat("Posterior distribution of K ...\n")
  Kd   <- posterior_K(out)
  if (verbose) cat(sprintf("  Mode K = %d,  Mean K = %.2f\n",
                           Kd$k[which.max(Kd$proportion)],
                           sum(Kd$k * Kd$proportion)))
  if (verbose) print(Kd)
  
  if (verbose) cat("Computing Posterior Similarity Matrix (may take a moment) ...\n")
  psm  <- compute_psm(out)
  
  if (verbose) cat("Binder loss partition ...\n")
  binder <- salso_partition(out, loss = loss)
  K_rep  <- length(unique(binder$partition))
  if (verbose) cat(sprintf("  Representative partition: K = %d clusters\n", K_rep))
  
  if (verbose) cat("Per-SNP BMA causal effects ...\n")
  bma  <- snp_bma_effects(out, ci)
  
  if (verbose) cat("Cluster-specific summaries ...\n")
  clus <- cluster_summaries(binder, bma, out, ci)
  if (verbose) print(clus)
  
  plei <- NULL
  if (!is.null(dat)) {
    if (verbose) cat("Pleiotropy check ...\n")
    plei <- pleiotropy_check(out, dat)
    if (verbose) cat(sprintf("  Posterior mean T = %.3f  (should be ~1 under correct spec)\n",
                             plei$mean_T))
    if (verbose) cat(sprintf("  95%% CI for T: [%.3f, %.3f]\n",
                             plei$ci_T[1], plei$ci_T[2]))
  }
  
  if (!is.null(fig_dir)) {
    dir.create(fig_dir, showWarnings = FALSE)
    if (verbose) cat("Producing figures ...\n")
    plot_psm(psm, binder, snp_names = snp_names,
             file = file.path(fig_dir, "fig_psm.png"),
             name_cex     = name_cex,
             width        = 1600,
             height       = 1500,
             res          = 220)  
    wald <- dat$Gamma_hat / dat$gamma_hat
    plot_snp_effects(bma, binder, snp_names = snp_names,
                     wald = wald, ordering = ordering,
                     file = file.path(fig_dir, "fig_snp_effects.png"),
                     name_cex = name_cex)
    plot_posterior_K(Kd, file = file.path(fig_dir, "fig_K_posterior.png"))
    plot_traces(out, pleiotropy = pleiotropy, file = file.path(fig_dir, "fig_traces.png"))
    if (verbose) cat("  Figures written to", fig_dir, "\n")
  }
  
  invisible(list(
    out        = out,
    diag       = diag,
    psm        = psm,
    binder     = binder,
    bma        = bma,
    clusters   = clus,
    K_posterior= Kd,
    pleiotropy = plei
  ))
}
