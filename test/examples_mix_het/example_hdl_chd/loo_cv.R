dat_post <- hdl_chd@data
dat_post[, 3] <- dat_post[, 3]^2
dat_post[, 4] <- dat_post[, 4]^2
names(dat_post) <- c("gamma_hat", "Gamma_hat", "sigma2X", "sigma2Y")
write.csv(dat_post, file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/hdl_chd.csv")

dir_MR <- "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR"
load(file.path(dir_MR, "out_MR.rda"))
write.table(results$out[["gamma"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/gamma_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
write.table(results$out[["psi"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/psi_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
write.table(rep(0, results$out[["S"]]),
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/tau_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
write.table(results$out[["xi"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/xi_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
beta_lst <- results$out[["beta"]]
max_K <- max(lengths(beta_lst))
beta_post <- t(sapply(beta_lst, function(x) {
  length(x) <- max_K
  x
}))
write.table(beta_post,
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/beta_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")

dir_PMR <- "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/PMR"
load(file.path(dir_PMR, "out_PMR.rda"))
write.table(results$out[["gamma"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/PMR/gamma_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
write.table(results$out[["psi"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/PMR/psi_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
write.table(results$out[["tau"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/PMR/tau_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
write.table(results$out[["xi"]],
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/PMR/xi_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")
beta_lst <- results$out[["beta"]]
max_K <- max(lengths(beta_lst))
beta_post <- t(sapply(beta_lst, function(x) {
  length(x) <- max_K
  x
}))
write.table(beta_post,
            file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/PMR/beta_post.csv",
            col.names = FALSE, row.names = FALSE, sep = ",")

# =============================================================================
# loo_compare.R  —  PSIS-LOO comparison: BNPM-MR vs BNPM-PMR
#
# Reads the posterior CSV files written by mr_mixture.cpp and computes
# Pareto-smoothed importance-sampling LOO-CV (PSIS-LOO) for each model,
# then produces a formal comparison with diagnostics and a WAIC cross-check.
#
# Reference: Vehtari, Gelman & Gabry (2017) "Practical Bayesian model
#   evaluation using leave-one-out cross-validation and WAIC."
#   Statistics and Computing 27(5): 1413–1432.
#
# Usage:
#   Rscript loo_compare.R  <data.csv>  <dir_MR>  <dir_PMR>  [burnin]
#
#   data.csv   CSV with columns: gamma_hat, Gamma_hat, sigma2X, sigma2Y
#              (same file passed to the C++ sampler)
#   dir_MR     directory containing BNPM-MR  posterior CSVs
#   dir_PMR    directory containing BNPM-PMR posterior CSVs
#   burnin     integer, number of leading draws to discard (default: 0)
#              Set this to your actual burn-in length. LOO on pre-convergence
#              draws is invalid.
#
# Expected files in each directory (no header, written by mr_mixture.cpp):
#   gamma_post.csv   S × 1
#   psi_post.csv     S × 1
#   tau_post.csv     S × 1   (all-zero for BNPM-MR)
#   xi_post.csv      S × p   comma-separated, 1-INDEXED cluster labels
#   beta_post.csv    S × maxK  comma-separated, NA-padded beyond K_s
# =============================================================================

suppressPackageStartupMessages({
  for (pkg in c("loo")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
})

# =============================================================================
# 1.  Log-likelihood — exact R translation of log_like() in mr_mixture.cpp
#
#     Computes the log joint density of (γ̂_j, Γ̂_j) after integrating out
#     the per-SNP confounder γ_j ~ N(γ, ψ²):
#
#       v    = ψ²(σ²_Y + τ² + β²σ²_X) + σ²_X(σ²_Y + τ²)
#       quad = (β²ψ² + σ²_Y + τ²)(γ̂ − γ)²
#              − 2βψ²(γ̂ − γ)(Γ̂ − γβ)
#              + (ψ² + σ²_X)(Γ̂ − γβ)²
#       ℓ    = −log(2π) − ½ log v − quad/(2v)
#
#     All arguments except gamma, psi2, tau2 may be length-p vectors.
#     Returns a length-p numeric vector.
# =============================================================================

ll_vec <- function(gh, Gh, s2X, s2Y, beta, gamma, psi2, tau2) {
  LOG2PI <- log(2 * pi)
  b2   <- beta^2
  syt2 <- s2Y + tau2                          # σ²_Y + τ²
  v    <- psi2 * (syt2 + b2 * s2X) + s2X * syt2
  dg   <- gh - gamma                          # γ̂ − γ
  dG   <- Gh - gamma * beta                   # Γ̂ − γβ
  quad <- (b2 * psi2 + syt2) * dg^2 -
           2 * beta * psi2    * dg * dG +
           (psi2 + s2X)       * dG^2
  out  <- -LOG2PI - 0.5 * log(v) - 0.5 * quad / v
  # Mirror C++ guard: if (v <= 0) return NEG_INF
  # In R, log(0)=-Inf and 0/0=NaN; replace both with -Inf explicitly.
  out[!is.finite(v) | v <= 0] <- -Inf
  out
}

# =============================================================================
# 2.  Unit test: verify ll_vec against a hand-computed reference value.
#     This would catch any transcription error vs the C++ code.
# =============================================================================

.validate_ll <- function() {
  LOG2PI <- log(2 * pi)
  
  # Arbitrary non-degenerate values
  gh <- 0.05; Gh <- 0.10; s2X <- 0.01; s2Y <- 0.04
  beta <- 0.5; gamma <- 0.02; psi2 <- 0.09; tau2 <- 0.0025

  # Compute by hand (same algebra as C++)
  b2   <- beta^2
  syt2 <- s2Y + tau2
  v    <- psi2 * (syt2 + b2 * s2X) + s2X * syt2
  dg   <- gh - gamma
  dG   <- Gh - gamma * beta
  quad <- (b2 * psi2 + syt2) * dg^2 -
           2 * beta * psi2    * dg * dG +
           (psi2 + s2X)       * dG^2
  expected <- -LOG2PI - 0.5 * log(v) - 0.5 * quad / v

  got <- ll_vec(gh, Gh, s2X, s2Y, beta, gamma, psi2, tau2)

  if (abs(got - expected) > 1e-12)
    stop(sprintf("ll_vec validation FAILED: got %.15f, expected %.15f", got, expected))

  # Also check that tau2=0 gives a different result from tau2>0 (sanity)
  ll0 <- ll_vec(gh, Gh, s2X, s2Y, beta, gamma, psi2, tau2 = 0)
  if (isTRUE(all.equal(ll0, got)))
    stop("ll_vec validation FAILED: tau2=0 and tau2>0 gave identical results")

  cat("ll_vec unit test: PASSED\n")
}
.validate_ll()

# =============================================================================
# 3.  Read posterior files
#
#     Returns a list with named elements:
#       gamma  numeric[S]       global confounding mean
#       psi    numeric[S]       psi (not squared)
#       tau    numeric[S]       tau (0 for every draw under BNPM-MR)
#       xi     integer[S × p]  cluster labels, 1-INDEXED (as written by C++)
#       beta   numeric[S × K]  cluster betas, NA-padded to maxK
#       S, p   integers
# =============================================================================

read_posterior <- function(dir, burnin = 0L) {

  read1 <- function(fname)
    scan(file.path(dir, fname), quiet = TRUE)

  gamma <- read1("gamma_post.csv")
  psi   <- read1("psi_post.csv")
  tau   <- read1("tau_post.csv")

  # xi: no header, entries are integers 1 … K_s
  xi <- as.matrix(read.csv(file.path(dir, "xi_post.csv"),
                            header = FALSE,
                            colClasses = "integer"))

  # beta: no header, NA-padded columns
  beta <- as.matrix(read.csv(file.path(dir, "beta_post.csv"),
                              header = FALSE, na.strings = "NA"))

  S_raw <- length(gamma)
  stopifnot(length(psi)  == S_raw,
            length(tau)  == S_raw,
            nrow(xi)     == S_raw,
            nrow(beta)   == S_raw)

  # Discard burn-in
  if (burnin >= S_raw)
    stop(sprintf("burnin (%d) >= total draws (%d)", burnin, S_raw))

  keep <- seq(burnin + 1L, S_raw)
  list(gamma = gamma[keep],
       psi   = psi[keep],
       tau   = tau[keep],
       xi    = xi[keep, , drop = FALSE],
       beta  = beta[keep, , drop = FALSE],
       S     = length(keep),
       p     = ncol(xi))
}

# =============================================================================
# 4.  Build log-likelihood matrix  ll[s, j]  (S × p)
#
#     ll[s, j] = log p(γ̂_j, Γ̂_j | θ^(s))
#
#     where θ^(s) = (γ^(s), β★^(s)_{ξ_j^(s)}, ψ^(s), τ^(s)).
#
#     NOTE ON CONDITIONAL LOO:
#     This uses ξ_j^(s) from the full-data posterior, i.e., it is the
#     *conditional* log-likelihood.  The PSIS-LOO approximation to the
#     true leave-one-out predictive density therefore conditions on the
#     cluster assignment of observation j.  This is mildly optimistic in
#     absolute terms.  However, because both BNPM-MR and BNPM-PMR share
#     the same DP mixture structure, the same approximation applies to
#     both models identically, so the *difference* in ELPD-LOO is an
#     unbiased comparison.  (See Vehtari et al. 2017, §5.)
# =============================================================================

build_ll_matrix <- function(post, data) {

  S <- post$S
  p <- post$p
  ll <- matrix(NA_real_, nrow = S, ncol = p)

  for (s in seq_len(S)) {

    # post$xi[s, ] is 1-indexed; use it to look up the beta for each SNP.
    # beta[s, xi[s,j]] is always a valid non-NA entry because xi[s,j] <= K_s.
    xi_s     <- post$xi[s, ]            # integer[p], values in 1 … K_s
    beta_j   <- post$beta[s, xi_s]      # numeric[p]

    if (anyNA(beta_j))
      stop(sprintf("NA in beta lookup at draw s=%d — xi/beta mismatch", s))

    ll[s, ] <- ll_vec(
      data$gamma_hat, data$Gamma_hat,
      data$sigma2X,   data$sigma2Y,
      beta_j,
      gamma = post$gamma[s],
      psi2  = post$psi[s]^2,
      tau2  = post$tau[s]^2
    )

    if (s %% 500 == 0L || s == S)
      cat(sprintf("  %5d / %d\r", s, S))
  }
  cat("\n")
  ll
}

# =============================================================================
# 5.  PSIS-LOO + WAIC for one model
# =============================================================================

compute_loo_waic <- function(ll, label) {

  n_neginf <- sum(!is.finite(ll))
  if (n_neginf > 0)
    warning(sprintf("[%s] %d non-finite ll entries — check model/data alignment",
                    label, n_neginf))

  S <- nrow(ll)

  # relative_eff: ESS-based MCMC correction for the IS variance estimator.
  # chain_id = rep(1, S) for a single chain.
  r_eff <- loo::relative_eff(exp(ll), chain_id = rep(1L, S))

  loo_obj  <- loo::loo(ll,  r_eff = r_eff)
  waic_obj <- loo::waic(ll)

  list(loo = loo_obj, waic = waic_obj)
}

# =============================================================================
# 6.  Main
# =============================================================================

data_path <- "hdl_chd.csv"
dir_mr    <- dir_MR
dir_pmr   <- dir_PMR
burnin    <- 0

# ── Load data ─────────────────────────────────────────────────────────────────

cat("Reading data:", data_path, "\n")
data <- read.csv("/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/hdl_chd.csv")
required <- c("gamma_hat", "Gamma_hat", "sigma2X", "sigma2Y")
missing  <- setdiff(required, names(data))
if (length(missing) > 0)
  stop("Missing columns in data CSV: ", paste(missing, collapse = ", "))

p <- nrow(data)
cat(sprintf("  p = %d SNPs\n", p))

# Sanity-check that sigma2X and sigma2Y are truly variances (not SEs).
# The C++ sampler expects VARIANCES.  If the typical value looks like a SE
# (e.g., median > 0.05), this warning fires.
if (median(data$sigma2X) > 0.05 || median(data$sigma2Y) > 0.05)
  warning("Median sigma2X or sigma2Y > 0.05 — are these variances or SEs?",
          " The sampler expects variances.")

# ── Load posteriors ────────────────────────────────────────────────────────────

cat("\nReading BNPM-MR posterior (burnin =", burnin, ") ...\n")
post_mr  <- read_posterior(dir_mr,  burnin)

cat("Reading BNPM-PMR posterior (burnin =", burnin, ") ...\n")
post_pmr <- read_posterior(dir_pmr, burnin)

stopifnot(post_mr$p == p, post_pmr$p == p)

cat(sprintf("  BNPM-MR  : %d post-burnin draws, max K = %d\n",
            post_mr$S,  ncol(post_mr$beta)))
cat(sprintf("  BNPM-PMR : %d post-burnin draws, max K = %d\n",
            post_pmr$S, ncol(post_pmr$beta)))

# Verify tau is truly zero under BNPM-MR
if (any(post_mr$tau != 0))
  stop("BNPM-MR tau_post contains non-zero values — wrong directory?")

# ── Build log-likelihood matrices ──────────────────────────────────────────────

cat("\nBuilding ll matrix: BNPM-MR ...\n")
ll_mr  <- build_ll_matrix(post_mr,  data)

cat("Building ll matrix: BNPM-PMR ...\n")
ll_pmr <- build_ll_matrix(post_pmr, data)

# ── PSIS-LOO and WAIC ─────────────────────────────────────────────────────────

cat("\nRunning PSIS-LOO and WAIC ...\n")
res_mr  <- compute_loo_waic(ll_mr,  "BNPM-MR")
res_pmr <- compute_loo_waic(ll_pmr, "BNPM-PMR")

# ── Report ─────────────────────────────────────────────────────────────────────

hr <- strrep("─", 60)

cat("\n", hr, "\n", sep = "")
cat("BNPM-MR  (pleiotropy = FALSE)\n")
cat(hr, "\n", sep = "")
print(res_mr$loo)

cat("\n", hr, "\n", sep = "")
cat("BNPM-PMR (pleiotropy = TRUE)\n")
cat(hr, "\n", sep = "")
print(res_pmr$loo)

cat("\n", hr, "\n", sep = "")
cat("PSIS-LOO model comparison\n")
cat("  Best model on first row; elpd_diff = 0 means it is the reference.\n")
cat("  A negative elpd_diff means that model is worse.\n")
cat("  Rule of thumb: |elpd_diff / se_diff| > 2 => meaningful difference.\n")
cat(hr, "\n", sep = "")
cmp <- loo_compare(res_mr$loo, res_pmr$loo)
print(cmp)

# Signed summary on elpd scale
elpd_mr  <- res_mr$loo$estimates["elpd_loo",  "Estimate"]
elpd_pmr <- res_pmr$loo$estimates["elpd_loo", "Estimate"]
se_mr    <- res_mr$loo$estimates["elpd_loo",  "SE"]
se_pmr   <- res_pmr$loo$estimates["elpd_loo", "SE"]

cat(sprintf("\n  BNPM-MR  ELPD-LOO : %+.2f  (SE %.2f)\n", elpd_mr,  se_mr))
cat(sprintf("  BNPM-PMR ELPD-LOO : %+.2f  (SE %.2f)\n", elpd_pmr, se_pmr))

# loo_compare rows come out ranked best-first; extract the named difference
diff_val <- cmp[, "elpd_diff"]
diff_se  <- cmp[, "se_diff"]
# The second row is always the worse model
better_name <- rownames(cmp)[1]
diff_magnitude <- abs(diff_val[2])
z_score <- diff_magnitude / diff_se[2]

cat(sprintf("\n  Better model (PSIS-LOO): %s\n", better_name))
cat(sprintf("  ELPD difference : %.2f  (SE %.2f,  z = %.2f)\n",
            diff_magnitude, diff_se[2], z_score))
if (z_score > 2) {
  cat("  => Difference is > 2 SE: considered MEANINGFUL.\n")
} else if (z_score > 1) {
  cat("  => Difference is 1–2 SE: WEAK evidence.\n")
} else {
  cat("  => Difference is < 1 SE: models are essentially EQUIVALENT by LOO.\n")
}

# ── Pareto k̂ diagnostics ──────────────────────────────────────────────────────

k_mr  <- res_mr$loo$diagnostics$pareto_k
k_pmr <- res_pmr$loo$diagnostics$pareto_k

cat("\n", hr, "\n", sep = "")
cat("Pareto k̂ diagnostics\n")
cat("  k̂ < 0.5  : good;  0.5–0.7 : OK;  > 0.7 : LOO estimate unreliable\n")
cat(hr, "\n", sep = "")
cat("  BNPM-MR  :", sum(k_mr  < 0.5), "good,",
                    sum(k_mr  >= 0.5 & k_mr  < 0.7), "OK,",
                    sum(k_mr  >= 0.7), "bad\n")
cat("  BNPM-PMR :", sum(k_pmr < 0.5), "good,",
                    sum(k_pmr >= 0.5 & k_pmr < 0.7), "OK,",
                    sum(k_pmr >= 0.7), "bad\n")

if (any(k_mr > 0.7))
  cat("  BNPM-MR  high-k̂ SNPs (1-indexed):",
      paste(which(k_mr > 0.7), collapse = ", "), "\n")
if (any(k_pmr > 0.7))
  cat("  BNPM-PMR high-k̂ SNPs (1-indexed):",
      paste(which(k_pmr > 0.7), collapse = ", "), "\n")

if (any(k_mr > 0.7) || any(k_pmr > 0.7)) {
  cat("\n  WARNING: High-k̂ observations mean the PSIS approximation is\n")
  cat("  unreliable for those SNPs. Consider moment matching (loo::loo_moment_match)\n")
  cat("  or exact LOO for the flagged SNPs to verify the comparison.\n")
}

# ── WAIC cross-check ──────────────────────────────────────────────────────────

cat("\n", hr, "\n", sep = "")
cat("WAIC (cross-check — should agree qualitatively with PSIS-LOO)\n")
cat(hr, "\n", sep = "")

waic_mr  <- res_mr$waic$estimates["elpd_waic",  "Estimate"]
waic_pmr <- res_pmr$waic$estimates["elpd_waic", "Estimate"]
se_waic_mr  <- res_mr$waic$estimates["elpd_waic",  "SE"]
se_waic_pmr <- res_pmr$waic$estimates["elpd_waic", "SE"]

cat(sprintf("  BNPM-MR  ELPD-WAIC : %+.2f  (SE %.2f)\n", waic_mr,  se_waic_mr))
cat(sprintf("  BNPM-PMR ELPD-WAIC : %+.2f  (SE %.2f)\n", waic_pmr, se_waic_pmr))

cmp_waic <- loo_compare(res_mr$waic, res_pmr$waic)
cat("WAIC comparison:\n")
print(cmp_waic)

if (sign(diff_val[2]) != sign(cmp_waic[2, "elpd_diff"])) {
  cat("\n  WARNING: PSIS-LOO and WAIC disagree on which model is better.\n")
} else {
  cat("\n  PSIS-LOO and WAIC agree on the ranking.\n")
}
