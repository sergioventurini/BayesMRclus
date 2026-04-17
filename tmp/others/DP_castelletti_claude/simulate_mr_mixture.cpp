// =============================================================================
// simulate_mr_mixture.cpp  —  Data simulator for the Bayesian MR Mixture Model
//
// Generates synthetic GWAS summary statistics under the mixture-of-causal-
// effects model assumed by bayesmr_mixture (mr_mixture.cpp).
//
// Compile:
//   g++ -O2 -std=c++17 -o simulate_mr_mixture simulate_mr_mixture.cpp
//
// Usage:
//   ./simulate_mr_mixture [options]
//
//   -p   <int>    number of SNPs              (default: 100)
//   -K   <int>    number of causal clusters   (default: 3)
//   -psi <float>  SD of true SNP-exposure effect γ_j   (default: 0.05)
//   -tau <float>  SD of directional pleiotropy u_j     (default: 0.00)
//   -n   <int>    effective GWAS sample size n (sets σ² ≈ 1/n)  (default: 50000)
//   -s   <int>    RNG seed                    (default: random)
//   -o   <str>    output prefix               (default: "sim")
//
// Output files:
//   <prefix>_data.csv       — input file for mr_mixture  (gamma_hat, Gamma_hat, sigma2X, sigma2Y)
//   <prefix>_truth.csv      — ground truth parameters    (gamma, beta, xi, u)
//
// =============================================================================
//
// Data-generating process
// -----------------------
//
//  1. Draw K cluster weights π from a symmetric Dirichlet(1/K).
//  2. Draw K true causal effects  β★_k ~ N(0, 1)  (well-separated clusters).
//  3. For each SNP j = 1 … p:
//       a. Assign cluster  ξ_j ~ Categorical(π).
//       b. Draw true SNP-exposure effect  γ_j ~ N(μ_γ, ψ²).
//       c. Draw pleiotropy noise  u_j ~ N(0, τ²).
//       d. Draw sampling variances  σ²_Xj, σ²_Yj ~ Uniform(0.5/n, 2.0/n).
//       e. Observe:
//            γ̂_j  | γ_j        ~  N(γ_j,                    σ²_Xj)
//            Γ̂_j  | γ_j, ξ_j  ~  N(γ_j · β★_{ξ_j} + u_j,  σ²_Yj)
//
// This matches the bivariate marginal in log_like() after integrating out γ_j.
// =============================================================================

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// =============================================================================
// 1.  Global RNG
// =============================================================================

static std::mt19937_64 rng;
inline void seed_rng(uint64_t s) { rng.seed(s); }

inline double rnorm_(double mu = 0.0, double sigma = 1.0) {
    return std::normal_distribution<double>(mu, sigma)(rng);
}
inline double runif_(double lo = 0.0, double hi = 1.0) {
    return std::uniform_real_distribution<double>(lo, hi)(rng);
}
inline double rgamma_(double shape, double rate) {
    return std::gamma_distribution<double>(shape, 1.0 / rate)(rng);
}
// Dirichlet(alpha, …, alpha) of length K
std::vector<double> rdirichlet_(int K, double concentration) {
    std::vector<double> x(K);
    double s = 0.0;
    for (int k = 0; k < K; ++k) { x[k] = rgamma_(concentration, 1.0); s += x[k]; }
    for (int k = 0; k < K; ++k) x[k] /= s;
    return x;
}
// Sample one index from a probability vector (already normalized)
int sample_discrete_(const std::vector<double>& w) {
    return std::discrete_distribution<int>(w.begin(), w.end())(rng);
}

// =============================================================================
// 2.  Command-line parsing helpers
// =============================================================================

struct Config {
    int         p           = 100;
    int         K           = 3;
    double      psi         = 0.05;
    double      tau         = 0.00;
    int         n           = 50000;
    uint64_t    seed        = 0;          // 0 → random
    std::string out_prefix  = "sim";
};

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "  -p  <int>    number of SNPs            (default: 100)\n"
        << "  -K  <int>    number of causal clusters (default: 3)\n"
        << "  -psi <float> SD of SNP-exposure effect (default: 0.05)\n"
        << "  -tau <float> SD of pleiotropy noise    (default: 0.00)\n"
        << "  -n  <int>    effective sample size     (default: 50000)\n"
        << "  -s  <int>    RNG seed (0 = random)     (default: random)\n"
        << "  -o  <str>    output prefix             (default: \"sim\")\n";
}

static Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];
        if ((key == "-p"   || key == "--p")   && i + 1 < argc) { cfg.p          = std::stoi(argv[++i]); }
        else if ((key == "-K"   || key == "--K")   && i + 1 < argc) { cfg.K          = std::stoi(argv[++i]); }
        else if ((key == "-psi" || key == "--psi") && i + 1 < argc) { cfg.psi        = std::stod(argv[++i]); }
        else if ((key == "-tau" || key == "--tau") && i + 1 < argc) { cfg.tau        = std::stod(argv[++i]); }
        else if ((key == "-n"   || key == "--n")   && i + 1 < argc) { cfg.n          = std::stoi(argv[++i]); }
        else if ((key == "-s"   || key == "--s")   && i + 1 < argc) { cfg.seed       = std::stoull(argv[++i]); }
        else if ((key == "-o"   || key == "--o")   && i + 1 < argc) { cfg.out_prefix = argv[++i]; }
        else if (key == "-h" || key == "--help") { usage(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown option: " << key << "\n"; usage(argv[0]); std::exit(1); }
    }
    return cfg;
}

// =============================================================================
// 3.  Simulator
// =============================================================================

int main(int argc, char* argv[]) {

    Config cfg = parse_args(argc, argv);

    // Seed
    const uint64_t seed = (cfg.seed == 0)
                        ? static_cast<uint64_t>(std::random_device{}())
                        : cfg.seed;
    seed_rng(seed);

    // ── Print simulation settings ─────────────────────────────────────────────
    std::cerr << "=== MR Mixture Data Simulator ===\n"
              << "  p   = " << cfg.p          << "  (SNPs)\n"
              << "  K   = " << cfg.K          << "  (clusters)\n"
              << "  psi = " << cfg.psi        << "  (SNP-exposure SD)\n"
              << "  tau = " << cfg.tau        << "  (pleiotropy SD)\n"
              << "  n   = " << cfg.n          << "  (effective sample size)\n"
              << "  seed= " << seed           << "\n"
              << "  out = " << cfg.out_prefix << "_data.csv / _truth.csv\n\n";

    const int    p   = cfg.p;
    const int    K   = cfg.K;
    const double psi = cfg.psi;
    const double tau = cfg.tau;
    const double n   = static_cast<double>(cfg.n);

    // ── Step 1: cluster mixing weights π ~ Dirichlet(1) ──────────────────────
    // Using concentration = 1 gives a uniform draw on the simplex, so cluster
    // sizes are variable and realistic (not all equal).
    const std::vector<double> pi_mix = rdirichlet_(K, 1.0);

    // ── Step 2: true causal effects β★_k ─────────────────────────────────────
    // We space K effects across a realistic range and then perturb slightly so
    // they are identifiable but not grid-like.
    // For K=1 the single effect is drawn N(0.5, 0.1²).
    std::vector<double> beta_star(K);
    if (K == 1) {
        beta_star[0] = rnorm_(0.5, 0.1);
    } else {
        // Evenly space K anchors in [-1, 1], then add small noise
        for (int k = 0; k < K; ++k) {
            const double anchor = -1.0 + 2.0 * k / (K - 1.0);
            beta_star[k] = anchor + rnorm_(0.0, 0.05);
        }
        // Shuffle so ordering is not trivially preserved
        std::shuffle(beta_star.begin(), beta_star.end(), rng);
    }

    // ── Step 3: simulate p SNPs ───────────────────────────────────────────────

    // Storage for the CSV columns
    std::vector<double> gamma_hat(p), Gamma_hat(p);
    std::vector<double> sigma2X(p),   sigma2Y(p);
    // Ground-truth columns
    std::vector<double> gamma_true(p), u_true(p);
    std::vector<int>    xi_true(p);
    std::vector<int>    cluster_sizes(K, 0);

    for (int j = 0; j < p; ++j) {

        // a) cluster assignment
        xi_true[j] = sample_discrete_(pi_mix);
        ++cluster_sizes[xi_true[j]];

        // b) true SNP-exposure effect  γ_j ~ N(0, ψ²)
        // We use μ_γ = 0 (the prior mean used in the sampler).
        gamma_true[j] = rnorm_(0.0, psi);

        // c) directional pleiotropy noise  u_j ~ N(0, τ²)
        u_true[j] = (tau > 0.0) ? rnorm_(0.0, tau) : 0.0;

        // d) sampling variances  σ²_Xj, σ²_Yj ~ Uniform(0.5/n, 2.0/n)
        // This mimics the spread seen in real GWAS where different SNPs are
        // tested in slightly different effective sample sizes.
        sigma2X[j] = runif_(0.5 / n, 2.0 / n);
        sigma2Y[j] = runif_(0.5 / n, 2.0 / n);

        // e) observed GWAS summary statistics
        //   γ̂_j  ~ N(γ_j,                       σ²_Xj)
        //   Γ̂_j  ~ N(γ_j · β★_{ξ_j} + u_j,     σ²_Yj)
        gamma_hat[j] = rnorm_(gamma_true[j],
                               std::sqrt(sigma2X[j]));
        Gamma_hat[j] = rnorm_(gamma_true[j] * beta_star[xi_true[j]] + u_true[j],
                               std::sqrt(sigma2Y[j]));
    }

    // ── Print summary ─────────────────────────────────────────────────────────
    std::cerr << "Cluster assignments:\n";
    for (int k = 0; k < K; ++k)
        std::cerr << "  cluster " << k + 1
                  << "  β★ = " << std::fixed << std::setprecision(4)
                  << beta_star[k]
                  << "  n_k = " << cluster_sizes[k]
                  << "  (π ≈ " << std::setprecision(3)
                  << pi_mix[k] << ")\n";
    std::cerr << "\n";

    // ── Write data CSV (input for mr_mixture) ─────────────────────────────────
    const std::string data_path  = cfg.out_prefix + "_data.csv";
    const std::string truth_path = cfg.out_prefix + "_truth.csv";

    {
        std::ofstream f(data_path);
        if (!f) { std::cerr << "Cannot open: " << data_path << "\n"; return 1; }
        f << "gamma_hat,Gamma_hat,sigma2X,sigma2Y\n";
        f << std::setprecision(10);
        for (int j = 0; j < p; ++j)
            f << gamma_hat[j] << ','
              << Gamma_hat[j] << ','
              << sigma2X[j]   << ','
              << sigma2Y[j]   << '\n';
    }

    // ── Write ground truth CSV ────────────────────────────────────────────────
    // Columns: snp, xi (1-indexed), beta_star, gamma, u
    {
        std::ofstream f(truth_path);
        if (!f) { std::cerr << "Cannot open: " << truth_path << "\n"; return 1; }
        f << "snp,xi,beta_star,gamma,u\n";
        f << std::setprecision(10);
        for (int j = 0; j < p; ++j)
            f << (j + 1)                           << ','   // 1-indexed
              << (xi_true[j] + 1)                  << ','   // 1-indexed
              << beta_star[xi_true[j]]             << ','
              << gamma_true[j]                     << ','
              << u_true[j]                         << '\n';
    }

    // ── Also write a small summary of the true β★ values ─────────────────────
    std::cerr << "True β★ values (1-indexed clusters):\n";
    for (int k = 0; k < K; ++k)
        std::cerr << "  β★[" << k + 1 << "] = " << std::fixed
                  << std::setprecision(6) << beta_star[k] << "\n";

    std::cerr << "\nWrote:\n"
              << "  " << data_path  << "  (" << p << " rows + header)\n"
              << "  " << truth_path << "  (" << p << " rows + header)\n";

    return 0;
}
