// =============================================================================
// diagnose_mr_mixture.cpp  —  Posterior diagnostics for bayesmr_mixture
//
// Reads the CSV files written by mr_mixture and the ground-truth CSV written
// by simulate_mr_mixture, then produces a comprehensive set of diagnostics
// that correctly account for the two fundamental difficulties of any DP
// mixture posterior:
//
//   (1) LABEL SWITCHING  — cluster labels are arbitrary and permuted across
//       iterations; it is meaningless to compare xi_post[s][j] directly with
//       xi_true[j], or beta_post[s][k] with the k-th true beta.
//
//   (2) VARIABLE K  — the number of occupied clusters changes across
//       iterations; scalar summaries of K must use its full posterior
//       distribution, not a single sample.
//
// Strategy
// --------
//  A. Label-switching-free per-SNP inference
//       beta_hat[s][j] = beta_post[s][ xi_post[s][j] - 1 ]
//     This quantity requires no cluster-label alignment.  Posterior mean and
//     95% credible intervals are compared directly to beta_true[j].
//
//  B. Posterior Similarity Matrix (PSM)
//       PSM[i][j] = (1/S_post) * Σ_s  1{ xi[s][i] == xi[s][j] }
//     The PSM is invariant to label switching and variable K.  It summarises
//     the entire clustering posterior in a single p×p matrix.
//
//  C. Point-estimate partition (medoid of PSM)
//     Among all post-burn-in samples, find the single partition ξ^(s*) whose
//     induced co-clustering matrix C^(s*) is closest (Frobenius distance) to
//     the PSM.  This is the least-squares representative sample under Binder
//     loss with equal misclassification costs.
//
//  D. Adjusted Rand Index (ARI)
//     Compares the medoid partition and the true partition without requiring
//     any label alignment.  ARI = 1 iff perfect recovery; ARI ≈ 0 for random.
//
//  E. Scalar parameters  (γ, ψ, τ, α, K)
//     Posterior mean ± SD, 95% equal-tail CI, and a coverage flag for the
//     true value.  K is treated as a discrete random variable.
//
// Compile:
//   g++ -O2 -std=c++17 -o diagnose_mr_mixture diagnose_mr_mixture.cpp
//
// Usage:
//   ./diagnose_mr_mixture <truth.csv> [burn_in_fraction] [posterior_prefix]
//
//   truth.csv          Output of simulate_mr_mixture: snp,xi,beta_star,gamma,u
//   burn_in_fraction   Fraction of iterations to discard (default: 0.5)
//   posterior_prefix   Path prefix for *_post.csv files (default: "")
//                      e.g. "results/" → reads results/gamma_post.csv etc.
//
// Output files (written to working directory):
//   diag_scalar.csv    Posterior summaries for γ, ψ, τ, α, K
//   diag_beta_snp.csv  Per-SNP: true β, posterior mean β̂, 95% CI, coverage
//   diag_psm.csv       p×p posterior similarity matrix (space-separated)
//   diag_clustering.csv  Per-SNP: true ξ, medoid ξ̂, correct flag
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <map>

// =============================================================================
// 1.  I/O helpers
// =============================================================================

// Read a single-column CSV (no header) of doubles.
std::vector<double> read_col(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::vector<double> v;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        v.push_back(std::stod(line));
    }
    return v;
}

// Read xi_post.csv  — rows of comma-separated integers (1-indexed).
// Returns S × p matrix stored row-major.  Converts to 0-indexed internally.
std::vector<std::vector<int>> read_xi(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::vector<std::vector<int>> mat;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        std::vector<int> row;
        while (std::getline(ss, tok, ','))
            row.push_back(std::stoi(tok) - 1);  // → 0-indexed
        mat.push_back(std::move(row));
    }
    return mat;
}

// Read beta_post.csv  — rows of comma-separated doubles or "NA".
// Returns S × K_s  (variable length per row).
std::vector<std::vector<double>> read_beta(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::vector<std::vector<double>> mat;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        std::vector<double> row;
        while (std::getline(ss, tok, ',')) {
            if (tok == "NA" || tok == "na") break;   // stop at first NA
            row.push_back(std::stod(tok));
        }
        mat.push_back(std::move(row));
    }
    return mat;
}

// Read truth CSV:  header  snp,xi,beta_star,gamma,u  (xi is 1-indexed)
struct TruthRow { int snp; int xi; double beta_star; double gamma; double u; };

std::vector<TruthRow> read_truth(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::string line;
    std::getline(f, line);  // skip header
    std::vector<TruthRow> rows;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        TruthRow r;
        std::getline(ss, tok, ','); r.snp       = std::stoi(tok);
        std::getline(ss, tok, ','); r.xi        = std::stoi(tok) - 1;  // 0-indexed
        std::getline(ss, tok, ','); r.beta_star = std::stod(tok);
        std::getline(ss, tok, ','); r.gamma     = std::stod(tok);
        std::getline(ss, tok, ','); r.u         = std::stod(tok);
        rows.push_back(r);
    }
    return rows;
}

// =============================================================================
// 2.  Summary-statistic helpers
// =============================================================================

struct Summary {
    double mean, sd, lo95, hi95;  // lo95/hi95 = 2.5th / 97.5th percentile
};

Summary summarise(std::vector<double> v) {
    const int n = static_cast<int>(v.size());
    assert(n > 1);
    double m = 0.0;
    for (double x : v) m += x;
    m /= n;
    double s2 = 0.0;
    for (double x : v) s2 += (x - m) * (x - m);
    s2 /= (n - 1);

    std::sort(v.begin(), v.end());
    const double lo = v[static_cast<int>(std::floor(0.025 * n))];
    const double hi = v[static_cast<int>(std::ceil (0.975 * n)) - 1];
    return {m, std::sqrt(s2), lo, hi};
}

// =============================================================================
// 3.  Adjusted Rand Index
// =============================================================================
//
//  Given two partitions a (length p) and b (length p) with arbitrary integer
//  labels, compute the ARI via the contingency-table formula.
//
//  ARI = [Σ_{ij} C(n_{ij},2) − t3] / [0.5*(t1+t2) − t3]
//  where t1 = Σ_i C(a_i,2),  t2 = Σ_j C(b_j,2),
//        t3 = t1*t2 / C(p,2),
//        n_{ij} = |{k: a[k]=i, b[k]=j}|.
//
double adjusted_rand_index(const std::vector<int>& a,
                           const std::vector<int>& b) {
    const int p = static_cast<int>(a.size());
    assert(p == static_cast<int>(b.size()));

    // Build contingency table via a map (handles arbitrary/non-contiguous labels)
    std::map<std::pair<int,int>, int> cont;
    std::map<int, int> ra, rb;
    for (int i = 0; i < p; ++i) {
        ++cont[{a[i], b[i]}];
        ++ra[a[i]];
        ++rb[b[i]];
    }

    auto c2 = [](long long n) -> double { return n <= 1 ? 0.0 : 0.5 * n * (n - 1); };

    double sum_nij = 0.0;
    for (auto& [key, n] : cont) sum_nij += c2(n);

    double t1 = 0.0;
    for (auto& [k, n] : ra) t1 += c2(n);
    double t2 = 0.0;
    for (auto& [k, n] : rb) t2 += c2(n);

    const double t3 = t1 * t2 / c2(p);
    const double denom = 0.5 * (t1 + t2) - t3;
    if (std::abs(denom) < 1e-12) return 1.0;  // both trivial or identical
    return (sum_nij - t3) / denom;
}

// =============================================================================
// 4.  Posterior Similarity Matrix (PSM)
// =============================================================================
//
//  Flat upper-triangle + diagonal stored row-major in a p×p dense matrix.
//  PSM[i*p + j] = proportion of post-burnin samples where xi[s][i] == xi[s][j].
//
std::vector<double> compute_psm(const std::vector<std::vector<int>>& xi_post,
                                int start, int p) {
    std::vector<double> psm(static_cast<size_t>(p) * p, 0.0);
    const int S_post = static_cast<int>(xi_post.size()) - start;
    for (int s = start; s < static_cast<int>(xi_post.size()); ++s) {
        const auto& xi = xi_post[s];
        for (int i = 0; i < p; ++i)
            for (int j = i; j < p; ++j)
                if (xi[i] == xi[j]) {
                    psm[i * p + j] += 1.0;
                    if (i != j) psm[j * p + i] += 1.0;
                }
    }
    const double inv = 1.0 / S_post;
    for (double& v : psm) v *= inv;
    return psm;
}

// =============================================================================
// 5.  Medoid partition (minimises Frobenius distance to PSM)
// =============================================================================
//
//  For each post-burnin sample s, compute:
//    D_s = Σ_{i≤j} (PSM[i][j] − 1{xi[s][i] == xi[s][j]})²
//  Return the ξ^(s*) with the smallest D_s.
//
//  This is the MAP-like representative of the clustering posterior under
//  Binder loss with equal misclassification costs, without any relabelling.
//
std::vector<int> medoid_partition(const std::vector<std::vector<int>>& xi_post,
                                  const std::vector<double>& psm,
                                  int start, int p) {
    double best_d = std::numeric_limits<double>::infinity();
    int    best_s = start;
    for (int s = start; s < static_cast<int>(xi_post.size()); ++s) {
        const auto& xi = xi_post[s];
        double d = 0.0;
        for (int i = 0; i < p; ++i) {
            for (int j = i; j < p; ++j) {
                const double diff = psm[i * p + j] - (xi[i] == xi[j] ? 1.0 : 0.0);
                d += diff * diff;
            }
        }
        if (d < best_d) { best_d = d; best_s = s; }
    }
    return xi_post[best_s];
}

// =============================================================================
// 6.  Main
// =============================================================================

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << " <truth.csv> [burn_in_fraction] [posterior_prefix]\n"
                  << "  truth.csv         from simulate_mr_mixture\n"
                  << "  burn_in_fraction  default 0.5\n"
                  << "  posterior_prefix  default \"\" (current dir)\n";
        return 1;
    }

    const std::string truth_path = argv[1];
    const double burn_frac = (argc >= 3) ? std::stod(argv[2]) : 0.5;
    const std::string pfx  = (argc >= 4) ? std::string(argv[3]) : std::string("");

    // ── Load ground truth ────────────────────────────────────────────────────
    const auto truth = read_truth(truth_path);
    if (truth.empty()) { std::cerr << "Empty truth file.\n"; return 1; }
    const int p = static_cast<int>(truth.size());

    // Per-SNP true causal effect (beta_star for its cluster)
    std::vector<double> beta_true(p);
    std::vector<int>    xi_true(p);
    for (int j = 0; j < p; ++j) {
        beta_true[j] = truth[j].beta_star;
        xi_true[j]   = truth[j].xi;   // 0-indexed
    }

    // True K
    const int K_true = 1 + *std::max_element(xi_true.begin(), xi_true.end());

    // True tau (if any u_j != 0 we infer tau was enabled; we can't recover
    // tau_true exactly from individual draws, so we use the empirical SD of u)
    double tau_true = 0.0;
    {
        double s2 = 0.0;
        for (const auto& r : truth) s2 += r.u * r.u;
        tau_true = std::sqrt(s2 / p);
    }

    std::cerr << "Truth loaded: p=" << p
              << "  K_true=" << K_true
              << "  tau_true(empirical)=" << std::fixed << std::setprecision(5)
              << tau_true << "\n";

    // ── Load posteriors ───────────────────────────────────────────────────────
    const auto gamma_post = read_col (pfx + "gamma_post.csv");
    const auto psi_post   = read_col (pfx + "psi_post.csv");
    const auto tau_post   = read_col (pfx + "tau_post.csv");
    const auto alpha_post = read_col (pfx + "alpha_post.csv");
    const auto xi_post    = read_xi  (pfx + "xi_post.csv");
    const auto beta_post  = read_beta(pfx + "beta_post.csv");

    if (gamma_post.empty() || xi_post.empty() || beta_post.empty()) {
        std::cerr << "Failed to load one or more posterior files.\n";
        return 1;
    }
    const int S = static_cast<int>(gamma_post.size());
    const int S_burn = static_cast<int>(std::ceil(burn_frac * S));
    const int S_post = S - S_burn;

    if (S_post <= 1) {
        std::cerr << "Too few post-burnin samples (S=" << S
                  << ", burn=" << S_burn << ").\n";
        return 1;
    }

    std::cerr << "Posterior: S=" << S << "  burn=" << S_burn
              << "  post-burnin=" << S_post << "\n\n";

    // ── Section A: Scalar parameters ─────────────────────────────────────────
    std::cout << "================================================================\n";
    std::cout << "SECTION A — Scalar parameter recovery\n";
    std::cout << "================================================================\n\n";

    // Helper: print a summary row
    auto print_scalar = [&](const std::string& name,
                            const std::vector<double>& full_chain,
                            double truth_val,
                            bool   truth_known = true) {
        std::vector<double> post(full_chain.begin() + S_burn, full_chain.end());
        const auto s = summarise(post);
        const bool covered = (s.lo95 <= truth_val && truth_val <= s.hi95);
        std::cout << std::left  << std::setw(8) << name
                  << "  mean=" << std::right << std::fixed << std::setprecision(5)
                  << std::setw(9) << s.mean
                  << "  sd=" << std::setw(7) << s.sd
                  << "  95%CI=[" << std::setw(9) << s.lo95
                  << ", " << std::setw(9) << s.hi95 << "]";
        if (truth_known)
            std::cout << "  truth=" << std::setw(9) << truth_val
                      << "  in_CI=" << (covered ? "YES" : "NO ");
        std::cout << "\n";
    };

    // γ — prior mean μ_γ = 0, so the truth here is 0 (the DGP uses μ_γ = 0)
    // More precisely, gamma_true in the DGP is the mean of SNP-level γ_j's;
    // the model parameter γ absorbs that mean, so a natural comparison is 0.
    print_scalar("gamma",   gamma_post,  0.0);
    print_scalar("psi",     psi_post,    0.0, false);   // truth not known from file
    print_scalar("tau",     tau_post,    tau_true);
    print_scalar("alpha",   alpha_post,  0.0, false);   // no point truth for alpha

    // K distribution
    {
        std::vector<double> K_chain(S_post);
        for (int s = 0; s < S_post; ++s)
            K_chain[s] = static_cast<double>(
                *std::max_element(xi_post[S_burn + s].begin(),
                                  xi_post[S_burn + s].end()) + 1);

        // Posterior distribution of K
        std::map<int, int> K_hist;
        for (double kd : K_chain) ++K_hist[static_cast<int>(kd)];

        const auto ks = summarise(K_chain);
        std::cout << std::left << std::setw(8) << "K"
                  << "  mean=" << std::fixed << std::setprecision(5)
                  << std::setw(9) << ks.mean
                  << "  sd=" << std::setw(7) << ks.sd
                  << "  95%CI=[" << std::setw(9) << ks.lo95
                  << ", " << std::setw(9) << ks.hi95 << "]"
                  << "  truth=" << std::setw(5) << K_true << "\n";

        std::cout << "\n  Posterior distribution of K:\n";
        for (auto& [k, cnt] : K_hist) {
            const double pct = 100.0 * cnt / S_post;
            const int bar = static_cast<int>(pct / 2.0);
            std::cout << "    K=" << std::setw(3) << k
                      << "  " << std::setw(6) << std::fixed << std::setprecision(2)
                      << pct << "%  |" << std::string(bar, '#') << "\n";
        }
    }

    // ── Section B: Per-SNP β recovery (label-switching free) ─────────────────
    std::cout << "\n================================================================\n";
    std::cout << "SECTION B — Per-SNP causal effect β recovery (label-switch free)\n";
    std::cout << "================================================================\n\n";
    std::cout << "  β̂_j^(s) = beta_post[s][ xi_post[s][j] ]   (no relabelling needed)\n\n";

    // For each SNP j, collect the per-sample causal effect
    std::vector<double> beta_pm(p, 0.0);      // posterior mean
    std::vector<double> beta_lo(p), beta_hi(p);
    std::vector<bool>   covered(p, false);

    for (int j = 0; j < p; ++j) {
        std::vector<double> draws(S_post);
        for (int s = 0; s < S_post; ++s) {
            const int ss = S_burn + s;
            const int clust = xi_post[ss][j];
            draws[s] = beta_post[ss][clust];
        }
        const auto sm = summarise(draws);
        beta_pm[j] = sm.mean;
        beta_lo[j] = sm.lo95;
        beta_hi[j] = sm.hi95;
        covered[j] = (sm.lo95 <= beta_true[j] && beta_true[j] <= sm.hi95);
    }

    // RMSE and coverage
    double rmse = 0.0;
    int    n_covered = 0;
    for (int j = 0; j < p; ++j) {
        rmse += (beta_pm[j] - beta_true[j]) * (beta_pm[j] - beta_true[j]);
        if (covered[j]) ++n_covered;
    }
    rmse = std::sqrt(rmse / p);
    const double coverage_pct = 100.0 * n_covered / p;

    std::cout << "  RMSE(β̂, β_true) = " << std::fixed << std::setprecision(6)
              << rmse << "\n";
    std::cout << "  95% CI coverage  = " << std::setprecision(1)
              << coverage_pct << "%  (" << n_covered << "/" << p
              << " SNPs)   [nominal: 95%]\n";

    // Print first 10 SNPs as a quick table
    std::cout << "\n  First " << std::min(p, 10) << " SNPs:\n";
    std::cout << "    " << std::setw(5) << "SNP"
              << std::setw(11) << "beta_true"
              << std::setw(11) << "post_mean"
              << std::setw(11) << "CI_lo"
              << std::setw(11) << "CI_hi"
              << std::setw(8)  << "covered"
              << "\n";
    for (int j = 0; j < std::min(p, 10); ++j)
        std::cout << "    " << std::setw(5) << (j + 1)
                  << std::fixed << std::setprecision(5)
                  << std::setw(11) << beta_true[j]
                  << std::setw(11) << beta_pm[j]
                  << std::setw(11) << beta_lo[j]
                  << std::setw(11) << beta_hi[j]
                  << std::setw(8)  << (covered[j] ? "YES" : "NO")
                  << "\n";

    // ── Section C: Posterior Similarity Matrix & clustering recovery ──────────
    std::cout << "\n================================================================\n";
    std::cout << "SECTION C — Clustering recovery via Posterior Similarity Matrix\n";
    std::cout << "================================================================\n\n";

    std::cerr << "Computing PSM (" << p << "x" << p
              << ", " << S_post << " post-burnin samples) ...\n";

    const auto psm = compute_psm(xi_post, S_burn, p);

    // PSM summary: average pairwise co-clustering for truly-same vs truly-different
    double psm_same = 0.0, psm_diff = 0.0;
    long   n_same   = 0,   n_diff   = 0;
    for (int i = 0; i < p; ++i) {
        for (int j = i + 1; j < p; ++j) {
            const double v = psm[i * p + j];
            if (xi_true[i] == xi_true[j]) { psm_same += v; ++n_same; }
            else                           { psm_diff += v; ++n_diff; }
        }
    }
    if (n_same > 0) psm_same /= n_same;
    if (n_diff > 0) psm_diff /= n_diff;

    std::cout << "  PSM summary (pairs sharing the same true cluster vs different):\n";
    std::cout << "    Mean PSM for truly co-clustered pairs  : "
              << std::fixed << std::setprecision(4) << psm_same
              << "  [ideal: 1.00]\n";
    std::cout << "    Mean PSM for truly separate-cluster pairs: "
              << std::setprecision(4) << psm_diff
              << "  [ideal: 0.00]\n";
    std::cout << "    Discrimination gap (same − diff)       : "
              << std::setprecision(4) << (psm_same - psm_diff) << "\n";

    // Medoid partition
    std::cerr << "Finding medoid partition ...\n";
    const auto xi_medoid = medoid_partition(xi_post, psm, S_burn, p);

    // Re-label medoid to 0-indexed compact integers
    {
        std::map<int, int> lmap;
        int next = 0;
        for (int xi : xi_medoid)
            if (lmap.find(xi) == lmap.end()) lmap[xi] = next++;
    }

    const int K_medoid = 1 + *std::max_element(xi_medoid.begin(), xi_medoid.end());
    const double ari = adjusted_rand_index(xi_medoid, xi_true);

    std::cout << "\n  Medoid partition (representative of the clustering posterior):\n";
    std::cout << "    K_medoid = " << K_medoid << "  (K_true = " << K_true << ")\n";
    std::cout << "    ARI(medoid, truth) = " << std::fixed << std::setprecision(4)
              << ari << "  [1.0 = perfect, ~0.0 = random]\n";

    // Cluster-level β★ summary for medoid partition
    // For each cluster in the medoid, report the posterior mean β̂ for its members
    // and the true β★ values of those members (should be homogeneous if K matches).
    {
        std::map<int, std::vector<int>> medoid_members;
        for (int j = 0; j < p; ++j) medoid_members[xi_medoid[j]].push_back(j);

        std::cout << "\n  Medoid clusters — estimated vs true β★ for each group:\n";
        std::cout << "    " << std::setw(10) << "cluster"
                  << std::setw(8)  << "size"
                  << std::setw(14) << "mean_beta_hat"
                  << std::setw(14) << "mean_beta_true"
                  << "\n";
        for (auto& [k, members] : medoid_members) {
            double mb_hat = 0.0, mb_true = 0.0;
            for (int j : members) { mb_hat += beta_pm[j]; mb_true += beta_true[j]; }
            mb_hat  /= members.size();
            mb_true /= members.size();
            std::cout << "    " << std::setw(10) << (k + 1)
                      << std::setw(8)  << members.size()
                      << std::fixed << std::setprecision(5)
                      << std::setw(14) << mb_hat
                      << std::setw(14) << mb_true
                      << "\n";
        }
    }

    // ── Section D: Binder loss trace (convergence check) ─────────────────────
    std::cout << "\n================================================================\n";
    std::cout << "SECTION D — Binder loss over iterations (convergence check)\n";
    std::cout << "================================================================\n\n";

    // Compute Binder loss of each post-burnin sample against the PSM.
    // A decreasing or stable trace is evidence of mixing; a declining trend
    // into the post-burnin period suggests burn-in was too short.
    const int n_report = std::min(S_post, 20);
    const int stride   = S_post / n_report;

    std::cout << "  Binder loss (lower = closer to posterior mean clustering):\n";
    double binder_min = std::numeric_limits<double>::infinity();
    double binder_max = 0.0;
    std::vector<double> binder_trace(S_post);
    for (int s = 0; s < S_post; ++s) {
        const auto& xi = xi_post[S_burn + s];
        double d = 0.0;
        for (int i = 0; i < p; ++i)
            for (int j = i + 1; j < p; ++j) {
                const double diff = psm[i * p + j] - (xi[i] == xi[j] ? 1.0 : 0.0);
                d += diff * diff;
            }
        binder_trace[s] = d;
        binder_min = std::min(binder_min, d);
        binder_max = std::max(binder_max, d);
    }
    for (int r = 0; r < n_report; ++r) {
        const int s = r * stride;
        const double pct = (binder_trace[s] - binder_min) /
                           (binder_max > binder_min ? binder_max - binder_min : 1.0);
        const int bar = static_cast<int>(40.0 * (1.0 - pct));
        std::cout << "    iter " << std::setw(6) << (S_burn + s + 1)
                  << "  loss=" << std::fixed << std::setprecision(2)
                  << std::setw(10) << binder_trace[s]
                  << "  |" << std::string(bar, '=')
                  << std::string(40 - bar, ' ') << "|\n";
    }

    // ── Write diagnostic CSVs ─────────────────────────────────────────────────
    std::cout << "\n================================================================\n";
    std::cout << "Writing output files ...\n";
    std::cout << "================================================================\n\n";

    // diag_scalar.csv
    {
        std::ofstream f("diag_scalar.csv");
        f << "param,mean,sd,ci_lo95,ci_hi95,truth,in_ci\n";
        auto row = [&](const std::string& name,
                       const std::vector<double>& chain,
                       double tv, bool has_truth) {
            std::vector<double> post(chain.begin() + S_burn, chain.end());
            const auto s = summarise(post);
            const bool cov = has_truth && (s.lo95 <= tv && tv <= s.hi95);
            f << name << "," << s.mean << "," << s.sd << ","
              << s.lo95 << "," << s.hi95 << ","
              << (has_truth ? std::to_string(tv) : "NA") << ","
              << (has_truth ? (cov ? "1" : "0") : "NA") << "\n";
        };
        row("gamma", gamma_post, 0.0, true);
        row("psi",   psi_post,   0.0, false);
        row("tau",   tau_post,   tau_true, true);
        row("alpha", alpha_post, 0.0, false);

        // K row
        std::vector<double> K_chain(S_post);
        for (int s = 0; s < S_post; ++s)
            K_chain[s] = static_cast<double>(
                *std::max_element(xi_post[S_burn+s].begin(),
                                  xi_post[S_burn+s].end()) + 1);
        const auto ks = summarise(K_chain);
        f << "K," << ks.mean << "," << ks.sd << ","
          << ks.lo95 << "," << ks.hi95 << "," << K_true << ","
          << (ks.lo95 <= K_true && K_true <= ks.hi95 ? "1" : "0") << "\n";
    }

    // diag_beta_snp.csv
    {
        std::ofstream f("diag_beta_snp.csv");
        f << "snp,xi_true,beta_true,beta_post_mean,ci_lo95,ci_hi95,in_ci\n";
        for (int j = 0; j < p; ++j)
            f << (j + 1) << ","
              << (xi_true[j] + 1) << ","
              << beta_true[j]  << ","
              << beta_pm[j]    << ","
              << beta_lo[j]    << ","
              << beta_hi[j]    << ","
              << (covered[j] ? 1 : 0) << "\n";
    }

    // diag_psm.csv  (space-separated, no header — can be read as matrix in R)
    {
        std::ofstream f("diag_psm.csv");
        f << std::fixed << std::setprecision(6);
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                f << psm[i * p + j];
                if (j + 1 < p) f << " ";
            }
            f << "\n";
        }
    }

    // diag_clustering.csv
    {
        std::ofstream f("diag_clustering.csv");
        f << "snp,xi_true,xi_medoid,correct\n";
        // "correct" here means medoid and truth agree on co-clustering with
        // all other SNPs, captured by ARI; per-SNP correctness still requires
        // a label match, so we use the contingency-table majority-vote mapping.
        // Build a best-label map from medoid → truth (greedy max-overlap).
        std::map<int, std::map<int,int>> cont;
        for (int j = 0; j < p; ++j) ++cont[xi_medoid[j]][xi_true[j]];
        std::map<int,int> best_map;
        for (auto& [km, inner] : cont) {
            int best_t = -1, best_n = -1;
            for (auto& [kt, n] : inner)
                if (n > best_n) { best_n = n; best_t = kt; }
            best_map[km] = best_t;
        }
        for (int j = 0; j < p; ++j)
            f << (j + 1) << ","
              << (xi_true[j] + 1) << ","
              << (xi_medoid[j] + 1) << ","
              << (best_map[xi_medoid[j]] == xi_true[j] ? 1 : 0) << "\n";
    }

    std::cout << "  diag_scalar.csv      — scalar parameter summaries\n";
    std::cout << "  diag_beta_snp.csv    — per-SNP β recovery\n";
    std::cout << "  diag_psm.csv         — " << p << "×" << p
              << " posterior similarity matrix\n";
    std::cout << "  diag_clustering.csv  — per-SNP true/medoid cluster labels\n";
    std::cout << "\n  ARI = " << std::fixed << std::setprecision(4) << ari
              << "  |  β RMSE = " << std::setprecision(6) << rmse
              << "  |  95% CI coverage = " << std::setprecision(1)
              << coverage_pct << "%\n\n";

    return 0;
}
