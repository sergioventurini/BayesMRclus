// =============================================================================
// mr_mixture.cpp  —  Bayesian Mendelian Randomisation Mixture Model
//
// Plain C++17 translation of the R implementation.
// No Rcpp, no Armadillo, no Eigen.
//
// Compile (recommended):
//   g++ -O3 -march=native -std=c++17 -o mr_mixture mr_mixture.cpp
//
// Usage:
//   ./mr_mixture <data.csv> <S> [seed]
//
//   data.csv  CSV file with header line:
//               gamma_hat,Gamma_hat,sigma2X,sigma2Y
//   S         number of MCMC iterations
//   seed      optional uint64 seed (default: non-deterministic)
//
// Output files (written to working directory):
//   gamma_post.csv   S × 1  – posterior draws of γ
//   psi_post.csv     S × 1  – posterior draws of ψ
//   tau_post.csv     S × 1  – posterior draws of τ (0 when pleiotropy=false)
//   alpha_post.csv   S × 1  – posterior draws of α (DP concentration)
//   xi_post.csv      S × p  – cluster indicators (1-indexed, comma-separated)
//   beta_post.csv    S × Kmax – β★ values (NA-padded for variable K)
//
// Algorithm:
//   Dirichlet Process mixture for heterogeneous causal effects, implemented
//   with Neal's Algorithm 8 (m auxiliary components) for the cluster update.
//   τ and ψ are updated via MH on the log scale; τ proposal is adaptive.
//   α is sampled exactly using the Escobar–West (1995) augmentation scheme.
// =============================================================================

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// =============================================================================
// 1.  Constants & global RNG
// =============================================================================

// log(2π)  — used in the bivariate Gaussian log-likelihood constant
static constexpr double LOG2PI  = 1.8378770664093454836;
static constexpr double NEG_INF = -std::numeric_limits<double>::infinity();

// A single 64-bit Mersenne Twister shared by all samplers.
// Seeded in main(); re-seed with seed_rng() if embedding this as a library.
static std::mt19937_64 rng;

inline void seed_rng(uint64_t s) { rng.seed(s); }

// =============================================================================
// 2.  Random-variate generators
// =============================================================================

inline double rnorm_(double mu = 0.0, double sigma = 1.0) {
    return std::normal_distribution<double>(mu, sigma)(rng);
}

inline double runif_() {
    return std::uniform_real_distribution<double>(0.0, 1.0)(rng);
}

// Gamma(shape, rate)  – note: std::gamma_distribution uses scale = 1/rate
inline double rgamma_(double shape, double rate) {
    return std::gamma_distribution<double>(shape, 1.0 / rate)(rng);
}

// Beta(a, b) via ratio of two Gamma draws
inline double rbeta_(double a, double b) {
    const double x = rgamma_(a, 1.0);
    return x / (x + rgamma_(b, 1.0));
}

// Uniform draw from {lo, lo+1, …, hi}  (inclusive)
inline int sample_int(int lo, int hi) {
    return std::uniform_int_distribution<int>(lo, hi)(rng);
}

// Sample one index from a (possibly unnormalized) probability vector.
// Uses std::discrete_distribution which internally normalizes.
inline int sample_discrete(const std::vector<double>& w) {
    return std::discrete_distribution<int>(w.begin(), w.end())(rng);
}

// =============================================================================
// 3.  Statistical helper functions
// =============================================================================

// ── 3a.  Log density of the half-t(ν, σ) distribution at x > 0 ──────────────
//
//   dht(x; ν, σ) = 2 · t_ν(x/σ) / σ
//
// where t_ν is the standard Student-t density with ν degrees of freedom.
// Corresponds to extraDistr::dht() in R.
double log_dht(double x, double nu, double sigma) {
    if (x <= 0.0) return NEG_INF;
    const double t = x / sigma;
    return std::log(2.0)
         + std::lgamma(0.5 * (nu + 1.0))
         - std::lgamma(0.5 * nu)
         - 0.5 * std::log(nu * M_PI)
         - 0.5 * (nu + 1.0) * std::log1p(t * t / nu)
         - std::log(sigma);
}

// ── 3b.  Log-likelihood for one observation (γ̂_j, Γ̂_j) ──────────────────────
//
// This is the bivariate-Gaussian marginal log-density obtained by integrating
// out γ_j.  The determinant of the 2×2 covariance matrix equals v.
//
//   v    = ψ²(σ²_Y + τ² + β²σ²_X) + σ²_X(σ²_Y + τ²)
//   quad = (β²ψ² + σ²_Y + τ²)(γ̂ − γ)²
//          − 2βψ²(γ̂ − γ)(Γ̂ − γβ)
//          + (ψ² + σ²_X)(Γ̂ − γβ)²
//
//   ℓ = −log(2π) − ½ log v − quad / (2v)
//
// Caller passes psi2 = ψ² and tau2 = τ² (pre-computed per iteration for speed).
inline double log_like(double gh,   double Gh,
                       double s2X,  double s2Y,
                       double beta, double gamma,
                       double psi2, double tau2) {
    const double b2    = beta * beta;
    const double syt2  = s2Y + tau2;          // σ²_Y + τ²  (re-used below)
    const double v     = psi2 * (syt2 + b2 * s2X) + s2X * syt2;
    if (v <= 0.0) return NEG_INF;
    const double dg    = gh - gamma;
    const double dG    = Gh - gamma * beta;
    const double quad  = (b2 * psi2 + syt2) * dg * dg
                       - 2.0 * beta * psi2    * dg * dG
                       + (psi2 + s2X)         * dG * dG;
    return -LOG2PI - 0.5 * std::log(v) - 0.5 * quad / v;
}

// ── 3c.  Log-sum-exp–stable normalization of log-weights → probabilities ──────
//
// Stores output in `out` (resized to match `lw`).
// Returns false if all entries are −∞ (all-zero output).
bool normalize_weights(const std::vector<double>& lw,
                       std::vector<double>&       out) {
    const int n = static_cast<int>(lw.size());
    out.resize(n);

    // Find the max of finite values (R code: max(prob[prob != -Inf]))
    double mx = NEG_INF;
    for (int i = 0; i < n; ++i)
        if (std::isfinite(lw[i]) && lw[i] > mx) mx = lw[i];

    if (!std::isfinite(mx)) {
        std::fill(out.begin(), out.end(), 0.0);
        return false;
    }

    double s = 0.0;
    for (int i = 0; i < n; ++i)
        if (std::isfinite(lw[i])) s += std::exp(lw[i] - mx);

    const double inv_s = 1.0 / s;
    for (int i = 0; i < n; ++i)
        out[i] = std::isfinite(lw[i]) ? std::exp(lw[i] - mx) * inv_s : 0.0;

    return true;
}

// ── 3d.  Conjugate Gaussian full conditional of γ ─────────────────────────────
//
// Prior:  γ ~ N(μ_γ, σ²_γ)
// Posterior precision A and mean B/A are accumulated by iterating over SNPs.
//
// In R the B numerator contains a telescoping cancellation:
//   β²ψ²γ̂ − β²ψ²γ̂ = 0
// leaving (σ²_Y + τ²)γ̂ + β σ²_X Γ̂  in the numerator.
double sample_gamma_fc(double mu_gamma,    double sigma2_gamma,
                       const std::vector<double>& gh,
                       const std::vector<double>& Gh,
                       const std::vector<double>& s2X,
                       const std::vector<double>& s2Y,
                       const std::vector<double>& beta_star,
                       const std::vector<int>&    xi,
                       double psi2, double tau2) {
    const int p = static_cast<int>(gh.size());
    double A = 1.0 / sigma2_gamma;
    double B = mu_gamma * A;

    for (int j = 0; j < p; ++j) {
        const double b    = beta_star[xi[j]];
        const double b2   = b * b;
        const double sx   = s2X[j];
        const double syt2 = s2Y[j] + tau2;       // σ²_Y + τ²
        const double v    = psi2 * (syt2 + b2 * sx) + sx * syt2;
        // Precision contribution: (σ²_Y + τ² + β²σ²_X) / v
        A += (syt2 + b2 * sx) / v;
        // Mean × precision contribution (after cancellation):
        //   [ (σ²_Y + τ²)γ̂ + β σ²_X Γ̂ ] / v
        B += (syt2 * gh[j] + b * sx * Gh[j]) / v;
    }
    return rnorm_(B / A, std::sqrt(1.0 / A));
}

// =============================================================================
// 4.  MCMC output container
// =============================================================================

struct MCMCOutput {
    std::vector<double>              gamma_post;   // S × 1
    std::vector<std::vector<double>> beta_post;    // S × K_s  (K varies)
    std::vector<double>              tau_post;     // S × 1  (0 if !pleiotropy)
    std::vector<double>              psi_post;     // S × 1
    std::vector<std::vector<int>>    xi_post;      // S × p  (0-indexed internally)
    std::vector<double>              alpha_post;   // S × 1
};

// =============================================================================
// 5.  Main MCMC sampler
// =============================================================================

// All vectors must have the same length p.
// xi_post entries are 0-indexed; convert to 1-indexed before writing to R/CSV.

MCMCOutput bayesmr_mixture(
    // ── Data ──────────────────────────────────────────────────────────────────
    const std::vector<double>& gamma_hat,
    const std::vector<double>& Gamma_hat,
    const std::vector<double>& sigma2X,
    const std::vector<double>& sigma2Y,
    // ── Hyperparameters ────────────────────────────────────────────────────────
    double mu_beta,    double sigma2_beta,
    double mu_gamma,   double sigma2_gamma,
    // ── Pleiotropy flag & half-t hyperparameters ───────────────────────────────
    bool   pleiotropy,
    double nu_tau,     double s_tau,   // ignored if pleiotropy == false
    double nu_psi,     double s_psi,
    // ── Proposal variances ─────────────────────────────────────────────────────
    double s2_betaprop,
    double s2_etaprop,    // adaptive if pleiotropy == true
    double s2_deltaprop,
    // ── Dirichlet Process ─────────────────────────────────────────────────────
    int    m,             // number of auxiliary components (Neal Alg 8)
    double a_alpha,       // Gamma(a_alpha, b_alpha) prior on DP concentration α
    double b_alpha,
    // ── MCMC ──────────────────────────────────────────────────────────────────
    int    S)
{
    const int p = static_cast<int>(gamma_hat.size());

    // ── 5a.  Initialisation ───────────────────────────────────────────────────

    int K = 2;
    std::vector<int> xi(p);

    // Sample xi ensuring all K clusters are occupied  (mirrors the R while loop)
    while (true) {
        for (int j = 0; j < p; ++j) xi[j] = sample_int(0, K - 1);
        bool c0 = false, c1 = false;
        for (int j = 0; j < p; ++j) { c0 |= (xi[j] == 0); c1 |= (xi[j] == 1); }
        if (c0 && c1) break;
    }

    std::vector<double> beta_star(K);
    for (int k = 0; k < K; ++k)
        beta_star[k] = rnorm_(mu_beta, std::sqrt(sigma2_beta));

    double tau   = 1e-2;
    double psi   = 1e-2;
    double alpha = 1.0;
    double gamma = 0.0;

    // ── 5b.  Allocate output storage ─────────────────────────────────────────

    MCMCOutput out;
    out.gamma_post.resize(S);
    out.beta_post.resize(S);
    out.tau_post.assign(S, 0.0);   // stays 0 if !pleiotropy
    out.psi_post.resize(S);
    out.xi_post.resize(S, std::vector<int>(p));
    out.alpha_post.resize(S);

    // Stores log(τ) draws; needed for the adaptive τ proposal
    std::vector<double> eta_post(S, std::log(tau));

    // ── 5c.  Pre-allocated working buffers (avoid repeated heap allocation) ───

    std::vector<double> log_w;          // log-weights for one SNP  (length K+m)
    std::vector<double> w_norm;         // normalised weights
    std::vector<double> beta_new;       // m auxiliary β draws
    std::vector<int>    xi_new(p);      // batch-updated cluster labels
    std::vector<int>    nk;             // cluster counts
    std::vector<bool>   used;           // which indices in [0, K+m) are occupied
    std::vector<int>    remap;          // old index → new (compacted) index
    std::vector<double> new_bstar;      // compacted β★ after relabelling

    // Cluster membership lists for the β★ MH update (avoids O(p) scan per k)
    std::vector<std::vector<int>> members;

    // ── 5d.  MCMC loop ────────────────────────────────────────────────────────

    std::cerr << "MCMC sampling ...\n";
    const int report_step = std::max(1, S / 100);

    for (int s = 0; s < S; ++s) {

        // Progress bar
        if ((s + 1) % report_step == 0 || s == S - 1)
            std::cerr << "\r  " << std::setw(3) << 100 * (s + 1) / S
                      << "%  [iter " << (s + 1) << "/" << S << "]   "
                      << std::flush;

        // Cache squared values — used extensively throughout this iteration
        const double psi2 = psi * psi;
        double       tau2 = tau * tau;   // mutable: updated after τ MH

        // ── Update γ  (conjugate Gibbs) ───────────────────────────────────────
        gamma = sample_gamma_fc(mu_gamma, sigma2_gamma,
                                gamma_hat, Gamma_hat,
                                sigma2X, sigma2Y,
                                beta_star, xi,
                                psi2, tau2);

        // ── Update ξ  (Neal's Algorithm 8, batch) ────────────────────────────
        K = static_cast<int>(beta_star.size());
        const int Km = K + m;

        // Draw m auxiliary β values from the DP base measure
        beta_new.resize(m);
        for (int k = 0; k < m; ++k)
            beta_new[k] = rnorm_(mu_beta, std::sqrt(sigma2_beta));

        // Cluster sizes from current xi (for leave-one-out counts below)
        nk.assign(K, 0);
        for (int j = 0; j < p; ++j) ++nk[xi[j]];

        log_w.resize(Km);
        const double log_a_m = std::log(alpha / static_cast<double>(m));

        for (int j = 0; j < p; ++j) {
            const int    xij = xi[j];
            const double ghj = gamma_hat[j], Ghj = Gamma_hat[j];
            const double sxj = sigma2X[j],   syj = sigma2Y[j];

            // Existing clusters (leave-one-out: subtract observation j)
            for (int k = 0; k < K; ++k) {
                const int n_jk = nk[k] - (xij == k ? 1 : 0);
                log_w[k] = (n_jk > 0)
                         ? std::log(static_cast<double>(n_jk)) +
                           log_like(ghj, Ghj, sxj, syj,
                                    beta_star[k], gamma, psi2, tau2)
                         : NEG_INF;
            }
            // Auxiliary new clusters
            for (int k = 0; k < m; ++k)
                log_w[K + k] = log_a_m +
                               log_like(ghj, Ghj, sxj, syj,
                                        beta_new[k], gamma, psi2, tau2);

            normalize_weights(log_w, w_norm);
            xi_new[j] = sample_discrete(w_norm);
        }
        xi = xi_new;  // batch update: all samples used the same old counts

        // ── Prune & relabel clusters ──────────────────────────────────────────
        // Identify which of the K + m candidate indices are actually occupied.
        used.assign(Km, false);
        for (int j = 0; j < p; ++j) used[xi[j]] = true;

        new_bstar.clear();
        remap.assign(Km, -1);
        for (int k = 0; k < Km; ++k) {
            if (!used[k]) continue;
            remap[k] = static_cast<int>(new_bstar.size());
            new_bstar.push_back(k < K ? beta_star[k] : beta_new[k - K]);
        }
        beta_star = std::move(new_bstar);
        K = static_cast<int>(beta_star.size());
        for (int j = 0; j < p; ++j) xi[j] = remap[xi[j]];

        // ── MH update for each β★_k (symmetric random walk) ──────────────────
        // Build per-cluster membership lists once; avoids O(p) scan per k.
        members.assign(K, {});
        for (int j = 0; j < p; ++j) members[xi[j]].push_back(j);

        const double sqrt_sbeta = std::sqrt(s2_betaprop);
        for (int k = 0; k < K; ++k) {
            const double bk   = beta_star[k];
            const double bkp  = rnorm_(bk, sqrt_sbeta);

            // Log prior ratio  (Normal prior on β★)
            double log_r = ((bk  - mu_beta) * (bk  - mu_beta)
                          - (bkp - mu_beta) * (bkp - mu_beta))
                         / (2.0 * sigma2_beta);

            // Log likelihood ratio over observations assigned to cluster k
            for (int j : members[k])
                log_r += log_like(gamma_hat[j], Gamma_hat[j],
                                  sigma2X[j],   sigma2Y[j],
                                  bkp, gamma, psi2, tau2)
                       - log_like(gamma_hat[j], Gamma_hat[j],
                                  sigma2X[j],   sigma2Y[j],
                                  bk,  gamma, psi2, tau2);

            if (std::log(runif_()) < log_r) beta_star[k] = bkp;
        }

        // ── MH update for τ = exp(η)  (if pleiotropy model) ──────────────────
        if (pleiotropy) {

            // Adaptive proposal: update s²_etaprop from empirical variance of η
            // after the first 10% of iterations  (mirrors R's adaptive step)
            if (s > static_cast<int>(0.1 * S) && s >= 2) {
                double me = 0.0;
                for (int i = 0; i < s; ++i) me += eta_post[i];
                me /= s;
                double ve = 0.0;
                for (int i = 0; i < s; ++i)
                    ve += (eta_post[i] - me) * (eta_post[i] - me);
                s2_etaprop = std::max(ve / (s - 1), 1e-8);
            }

            const double eta      = std::log(tau);
            const double eta_prop = rnorm_(eta, std::sqrt(s2_etaprop));
            const double tau_prop = std::exp(eta_prop);
            const double tau2_p   = tau_prop * tau_prop;

            // Log acceptance ratio:
            //   log π(τ') + η'  −  log π(τ) − η       (Jacobian for log-scale)
            //   + Σ_j [ ℓ(τ') − ℓ(τ) ]
            double log_r = log_dht(tau_prop, nu_tau, s_tau) + eta_prop
                         - log_dht(tau,      nu_tau, s_tau) - eta;
            for (int j = 0; j < p; ++j)
                log_r += log_like(gamma_hat[j], Gamma_hat[j],
                                  sigma2X[j],   sigma2Y[j],
                                  beta_star[xi[j]], gamma, psi2, tau2_p)
                       - log_like(gamma_hat[j], Gamma_hat[j],
                                  sigma2X[j],   sigma2Y[j],
                                  beta_star[xi[j]], gamma, psi2, tau2);

            if (std::log(runif_()) < log_r) {
                tau  = tau_prop;
                tau2 = tau2_p;   // keep tau2 in sync for the ψ update below
            }
            eta_post[s]     = std::log(tau);
            out.tau_post[s] = tau;
        }

        // ── MH update for ψ = exp(δ) ──────────────────────────────────────────
        {
            const double delta      = std::log(psi);
            const double delta_prop = rnorm_(delta, std::sqrt(s2_deltaprop));
            const double psi_prop   = std::exp(delta_prop);
            const double psi2_p     = psi_prop * psi_prop;

            double log_r = log_dht(psi_prop, nu_psi, s_psi) + delta_prop
                         - log_dht(psi,      nu_psi, s_psi) - delta;
            for (int j = 0; j < p; ++j)
                log_r += log_like(gamma_hat[j], Gamma_hat[j],
                                  sigma2X[j],   sigma2Y[j],
                                  beta_star[xi[j]], gamma, psi2_p, tau2)
                       - log_like(gamma_hat[j], Gamma_hat[j],
                                  sigma2X[j],   sigma2Y[j],
                                  beta_star[xi[j]], gamma, psi2,   tau2);

            if (std::log(runif_()) < log_r) psi = psi_prop;
        }

        // ── Update α  (Escobar & West 1995 exact Gamma mixture sampler) ───────
        //
        // Draw η_dp ~ Beta(α+1, p), then α from a mixture of two Gammas:
        //   α | K, η_dp  ∝  w₁ · Γ(a+K,   b−log η_dp)
        //                 + w₂ · Γ(a+K−1, b−log η_dp)
        // with mixing weights w₁ ∝ (a+K−1),  w₂ ∝ p(b−log η_dp).
        {
            const double eta_dp = rbeta_(alpha + 1.0, static_cast<double>(p));
            const double rate   = b_alpha - std::log(eta_dp);   // > 0 a.s.
            const double a1     = rgamma_(a_alpha + K,       rate);
            const double a2     = rgamma_(a_alpha + K - 1.0, rate);
            const double w1     = a_alpha + K - 1.0;
            const double w2     = static_cast<double>(p) * rate;
            alpha = (runif_() < w1 / (w1 + w2)) ? a1 : a2;
        }

        // ── Store samples ─────────────────────────────────────────────────────
        out.xi_post[s]    = xi;
        out.gamma_post[s] = gamma;
        out.beta_post[s]  = beta_star;
        out.alpha_post[s] = alpha;
        out.psi_post[s]   = psi;
    }

    std::cerr << "\nDone.\n";
    return out;
}

// =============================================================================
// 6.  I/O helpers
// =============================================================================

// Strip a trailing \r from a string in-place (handles Windows CRLF files).
static void strip_cr(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

// Split a (already \r-stripped) line on commas.
static std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream ss(line);
    std::string tok;
    while (std::getline(ss, tok, ',')) { strip_cr(tok); tokens.push_back(tok); }
    return tokens;
}

// Read a CSV file whose header contains (in any order) the four required
// column names: gamma_hat, Gamma_hat, sigma2X, sigma2Y.
// Scientific notation (e.g. 9.55E-05) and Windows line endings are handled.
bool read_csv(const std::string& path,
              std::vector<double>& gamma_hat,
              std::vector<double>& Gamma_hat,
              std::vector<double>& sigma2X,
              std::vector<double>& sigma2Y) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return false; }

    // ── Parse header → find column indices ────────────────────────────────
    std::string line;
    if (!std::getline(f, line)) {
        std::cerr << "Empty file: " << path << "\n"; return false;
    }
    strip_cr(line);
    const auto header = split_csv(line);

    int col_gh = -1, col_Gh = -1, col_sX = -1, col_sY = -1;
    for (int c = 0; c < static_cast<int>(header.size()); ++c) {
        const std::string& h = header[c];
        if (h == "gamma_hat") col_gh = c;
        else if (h == "Gamma_hat") col_Gh = c;
        else if (h == "sigma2X")   col_sX = c;
        else if (h == "sigma2Y")   col_sY = c;
    }
    if (col_gh < 0 || col_Gh < 0 || col_sX < 0 || col_sY < 0) {
        std::cerr << "Header must contain: gamma_hat, Gamma_hat, sigma2X, sigma2Y\n"
                  << "Found: ";
        for (const auto& h : header) std::cerr << h << " ";
        std::cerr << "\n";
        return false;
    }
    const int ncol_needed = std::max({col_gh, col_Gh, col_sX, col_sY}) + 1;

    // ── Read data rows ─────────────────────────────────────────────────────
    while (std::getline(f, line)) {
        strip_cr(line);
        if (line.empty()) continue;
        const auto toks = split_csv(line);
        if (static_cast<int>(toks.size()) < ncol_needed) {
            std::cerr << "Malformed row (need " << ncol_needed
                      << " columns): " << line << "\n";
            return false;
        }
        gamma_hat.push_back(std::stod(toks[col_gh]));
        Gamma_hat.push_back(std::stod(toks[col_Gh]));
        sigma2X  .push_back(std::stod(toks[col_sX]));
        sigma2Y  .push_back(std::stod(toks[col_sY]));
    }
    return !gamma_hat.empty();
}

// Write a 1-D vector of doubles, one per line
void write_vec(const std::string& path, const std::vector<double>& v) {
    std::ofstream f(path);
    f << std::setprecision(10);
    for (double x : v) f << x << "\n";
}

// Write xi_post (0-indexed internally → 1-indexed in output for R compatibility)
void write_xi(const std::string& path,
              const std::vector<std::vector<int>>& xi_post) {
    std::ofstream f(path);
    for (const auto& row : xi_post) {
        for (int j = 0; j < static_cast<int>(row.size()); ++j) {
            f << row[j] + 1;   // 1-indexed
            if (j + 1 < static_cast<int>(row.size())) f << ',';
        }
        f << '\n';
    }
}

// Write beta_post: variable-length rows, NA-padded to max K across all iters
void write_beta(const std::string& path,
                const std::vector<std::vector<double>>& beta_post) {
    int maxK = 0;
    for (const auto& b : beta_post)
        maxK = std::max(maxK, static_cast<int>(b.size()));

    std::ofstream f(path);
    f << std::setprecision(10);
    for (const auto& b : beta_post) {
        for (int k = 0; k < maxK; ++k) {
            if (k < static_cast<int>(b.size())) f << b[k];
            else                                f << "NA";
            if (k + 1 < maxK) f << ',';
        }
        f << '\n';
    }
}

// =============================================================================
// 7.  Entry point
// =============================================================================

int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <data.csv> <S> [seed]\n"
                  << "  data.csv  CSV with header containing (any order):\n"
                  << "              gamma_hat, Gamma_hat, sigma2X, sigma2Y\n"
                  << "  S         number of MCMC iterations\n"
                  << "  seed      optional RNG seed (default: random)\n";
        return 1;
    }

    const std::string data_path = argv[1];
    const int         S         = std::stoi(argv[2]);
    const uint64_t    seed      = (argc >= 4)
                                ? std::stoull(argv[3])
                                : static_cast<uint64_t>(std::random_device{}());
    seed_rng(seed);
    std::cerr << "Seed: " << seed << "\n";

    // Load data
    std::vector<double> gamma_hat, Gamma_hat, sigma2X, sigma2Y;
    if (!read_csv(data_path, gamma_hat, Gamma_hat, sigma2X, sigma2Y)) return 1;
    std::cerr << "Loaded p = " << gamma_hat.size() << " SNPs.\n";

    // ── Hyperparameters (mirrors R function defaults) ─────────────────────────
    const double mu_beta      = 0.0,   sigma2_beta  = 10.0;
    const double mu_gamma     = 0.0,   sigma2_gamma = 10.0;
    const bool   pleiotropy   = true; // set to true to enable τ updates
    const double nu_tau       = 3.0,   s_tau        = 0.01;  // used iff pleiotropy
    const double nu_psi       = 2.0,   s_psi        = 0.01;
    const double s2_betaprop  = 0.2,   s2_etaprop   = 0.2,  s2_deltaprop = 0.2;
    const int    m            = 2;
    const double a_alpha      = 3.0,   b_alpha      = 1.0;

    const auto t0 = std::chrono::high_resolution_clock::now();

    MCMCOutput out = bayesmr_mixture(
        gamma_hat, Gamma_hat, sigma2X, sigma2Y,
        mu_beta,  sigma2_beta,
        mu_gamma, sigma2_gamma,
        pleiotropy, nu_tau, s_tau, nu_psi, s_psi,
        s2_betaprop, s2_etaprop, s2_deltaprop,
        m, a_alpha, b_alpha, S);

    const auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();
    std::cerr << std::fixed << std::setprecision(3)
              << "Elapsed: " << elapsed << " s  ("
              << elapsed / S * 1000.0 << " ms / iteration)\n";

    // Write posterior samples
    write_vec ("gamma_post.csv",  out.gamma_post);
    write_vec ("alpha_post.csv",  out.alpha_post);
    write_vec ("psi_post.csv",    out.psi_post);
    write_vec ("tau_post.csv",    out.tau_post);
    write_xi  ("xi_post.csv",     out.xi_post);
    write_beta("beta_post.csv",   out.beta_post);

    std::cerr << "Posterior samples written to *_post.csv\n";
    return 0;
}
