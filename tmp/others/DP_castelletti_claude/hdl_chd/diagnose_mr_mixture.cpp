// =============================================================================
// diagnose_mr_mixture.cpp  —  Posterior diagnostics for bayesmr_mixture
//
// Works in two modes:
//
//   SIMULATION mode  (truth file provided)
//     Compares the posterior to ground truth from simulate_mr_mixture.
//     Computes β RMSE, 95%-CI coverage, ARI of medoid vs true partition.
//
//   REAL DATA mode   (no truth file, or truth = "-")
//     Reports only quantities derivable from the posterior alone:
//     scalar summaries (γ, ψ, τ, α, K), per-SNP β credible intervals,
//     Posterior Similarity Matrix, medoid partition, Binder loss trace.
//
// Compile:
//   g++ -O2 -std=c++17 -o diagnose_mr_mixture diagnose_mr_mixture.cpp
//
// Usage:
//   ./diagnose_mr_mixture [options]
//
//   --truth  <file>   Ground truth CSV from simulate_mr_mixture (optional;
//                     omit or pass "-" for real data with no ground truth)
//   --burn   <float>  Burn-in fraction (default 0.5)
//   --prefix <str>    Path prefix for *_post.csv files (default "")
//   --out    <str>    Output file prefix (default "diag")
//
// Output files:
//   <out>_scalar.csv      Posterior summaries for gamma, psi, tau, alpha, K
//   <out>_beta_snp.csv    Per-SNP posterior mean beta and 95% CI
//   <out>_psm.csv         p x p posterior similarity matrix (space-separated)
//   <out>_clustering.csv  Per-SNP medoid cluster label (+ truth if available)
//   <out>_binder.csv      Binder loss at every post-burnin iteration
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

// =============================================================================
// 1.  I/O helpers
// =============================================================================

static void strip_cr(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

std::vector<double> read_col(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::vector<double> v;
    std::string line;
    while (std::getline(f, line)) {
        strip_cr(line);
        if (line.empty()) continue;
        v.push_back(std::stod(line));
    }
    return v;
}

// xi_post.csv: rows of comma-separated 1-indexed integers -> 0-indexed
std::vector<std::vector<int>> read_xi(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::vector<std::vector<int>> mat;
    std::string line;
    while (std::getline(f, line)) {
        strip_cr(line);
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        std::vector<int> row;
        while (std::getline(ss, tok, ',')) {
            strip_cr(tok);
            row.push_back(std::stoi(tok) - 1);
        }
        mat.push_back(std::move(row));
    }
    return mat;
}

// beta_post.csv: variable-length rows, stop at first "NA"
std::vector<std::vector<double>> read_beta(const std::string& path) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open: " << path << "\n"; return {}; }
    std::vector<std::vector<double>> mat;
    std::string line;
    while (std::getline(f, line)) {
        strip_cr(line);
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        std::vector<double> row;
        while (std::getline(ss, tok, ',')) {
            strip_cr(tok);
            if (tok == "NA" || tok == "na") break;
            row.push_back(std::stod(tok));
        }
        mat.push_back(std::move(row));
    }
    return mat;
}

struct TruthRow { int snp, xi; double beta_star, gamma, u; };

// Returns empty vector when path == "-" or file is missing.
std::vector<TruthRow> read_truth(const std::string& path) {
    if (path == "-") return {};
    std::ifstream f(path);
    if (!f) {
        std::cerr << "Note: cannot open truth file '" << path
                  << "' — running in real-data mode.\n";
        return {};
    }
    std::string line;
    std::getline(f, line);   // skip header
    std::vector<TruthRow> rows;
    while (std::getline(f, line)) {
        strip_cr(line);
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        TruthRow r;
        std::getline(ss, tok, ','); r.snp       = std::stoi(tok);
        std::getline(ss, tok, ','); r.xi        = std::stoi(tok) - 1; // 0-indexed
        std::getline(ss, tok, ','); r.beta_star = std::stod(tok);
        std::getline(ss, tok, ','); r.gamma     = std::stod(tok);
        std::getline(ss, tok, ','); r.u         = std::stod(tok);
        rows.push_back(r);
    }
    return rows;
}

// =============================================================================
// 2.  Statistical helpers
// =============================================================================

struct Summary { double mean, sd, lo95, hi95; };

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
    const double hi = v[std::min(n - 1, static_cast<int>(std::ceil(0.975 * n)) - 1)];
    return {m, std::sqrt(s2), lo, hi};
}

// =============================================================================
// 3.  Adjusted Rand Index
// =============================================================================

double adjusted_rand_index(const std::vector<int>& a,
                           const std::vector<int>& b) {
    const int p = static_cast<int>(a.size());
    std::map<std::pair<int,int>, int> cont;
    std::map<int, int> ra, rb;
    for (int i = 0; i < p; ++i) {
        ++cont[{a[i], b[i]}]; ++ra[a[i]]; ++rb[b[i]];
    }
    auto c2 = [](long long n) -> double { return n <= 1 ? 0.0 : 0.5*n*(n-1); };
    double sij = 0.0;
    for (auto& [k, n] : cont) sij += c2(n);
    double t1 = 0.0; for (auto& [k, n] : ra) t1 += c2(n);
    double t2 = 0.0; for (auto& [k, n] : rb) t2 += c2(n);
    const double t3 = t1 * t2 / c2(p);
    const double den = 0.5*(t1+t2) - t3;
    if (std::abs(den) < 1e-12) return 1.0;
    return (sij - t3) / den;
}

// =============================================================================
// 4.  Posterior Similarity Matrix
// =============================================================================

std::vector<double> compute_psm(const std::vector<std::vector<int>>& xi_post,
                                int start, int p) {
    std::vector<double> psm(static_cast<size_t>(p) * p, 0.0);
    const int S_post = static_cast<int>(xi_post.size()) - start;
    for (int s = start; s < static_cast<int>(xi_post.size()); ++s) {
        const auto& xi = xi_post[s];
        for (int i = 0; i < p; ++i)
            for (int j = i; j < p; ++j)
                if (xi[i] == xi[j]) {
                    psm[i*p+j] += 1.0;
                    if (i != j) psm[j*p+i] += 1.0;
                }
    }
    const double inv = 1.0 / S_post;
    for (double& v : psm) v *= inv;
    return psm;
}

// =============================================================================
// 5.  Medoid partition
// =============================================================================

std::vector<int> medoid_partition(const std::vector<std::vector<int>>& xi_post,
                                  const std::vector<double>& psm,
                                  int start, int p) {
    double best = std::numeric_limits<double>::infinity();
    int    best_s = start;
    for (int s = start; s < static_cast<int>(xi_post.size()); ++s) {
        const auto& xi = xi_post[s];
        double d = 0.0;
        for (int i = 0; i < p; ++i)
            for (int j = i; j < p; ++j) {
                double diff = psm[i*p+j] - (xi[i] == xi[j] ? 1.0 : 0.0);
                d += diff * diff;
            }
        if (d < best) { best = d; best_s = s; }
    }
    return xi_post[best_s];
}

// =============================================================================
// 6.  Command-line parsing
// =============================================================================

struct Config {
    std::string truth_path  = "-";
    double      burn_frac   = 0.5;
    std::string post_prefix = "";
    std::string out_prefix  = "diag";
};

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " [options]\n"
        << "  --truth  <file>   truth CSV from simulate_mr_mixture\n"
        << "                    (omit or \"-\" for real data)\n"
        << "  --burn   <float>  burn-in fraction       (default 0.5)\n"
        << "  --prefix <str>    prefix for *_post.csv  (default \"\")\n"
        << "  --out    <str>    output file prefix      (default \"diag\")\n";
}

static Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string k = argv[i];
        if      (k == "--truth"  && i+1 < argc) cfg.truth_path  = argv[++i];
        else if (k == "--burn"   && i+1 < argc) cfg.burn_frac   = std::stod(argv[++i]);
        else if (k == "--prefix" && i+1 < argc) cfg.post_prefix = argv[++i];
        else if (k == "--out"    && i+1 < argc) cfg.out_prefix  = argv[++i];
        else if (k == "-h" || k == "--help") { usage(argv[0]); std::exit(0); }
        else { std::cerr << "Unknown option: " << k << "\n"; usage(argv[0]); std::exit(1); }
    }
    return cfg;
}

// =============================================================================
// 7.  Main
// =============================================================================

int main(int argc, char* argv[]) {

    const Config cfg  = parse_args(argc, argv);
    const auto&  pfx  = cfg.post_prefix;
    const auto&  opfx = cfg.out_prefix;

    // ── Load posteriors ───────────────────────────────────────────────────────
    const auto gamma_post = read_col (pfx + "gamma_post.csv");
    const auto psi_post   = read_col (pfx + "psi_post.csv");
    const auto tau_post   = read_col (pfx + "tau_post.csv");
    const auto alpha_post = read_col (pfx + "alpha_post.csv");
    const auto xi_post    = read_xi  (pfx + "xi_post.csv");
    const auto beta_post  = read_beta(pfx + "beta_post.csv");

    if (gamma_post.empty() || xi_post.empty() || beta_post.empty()) {
        std::cerr << "Failed to load posterior files from: "
                  << (pfx.empty() ? "current directory" : pfx) << "\n";
        return 1;
    }

    const int S      = static_cast<int>(gamma_post.size());
    const int p      = static_cast<int>(xi_post[0].size());
    const int S_burn = static_cast<int>(std::ceil(cfg.burn_frac * S));
    const int S_post = S - S_burn;

    if (S_post <= 1) {
        std::cerr << "Too few post-burnin samples (S=" << S
                  << ", burn=" << S_burn << ").\n"; return 1;
    }

    // ── Load truth (optional) ─────────────────────────────────────────────────
    const auto truth    = read_truth(cfg.truth_path);
    const bool truth_ok = !truth.empty() &&
                          static_cast<int>(truth.size()) == p;

    if (!truth.empty() && !truth_ok)
        std::cerr << "Warning: truth file has " << truth.size()
                  << " rows but posterior has p=" << p
                  << " — skipping truth comparisons.\n\n";

    std::vector<double> beta_true;
    std::vector<int>    xi_true;
    int    K_true   = -1;
    double tau_true = 0.0;

    if (truth_ok) {
        beta_true.resize(p); xi_true.resize(p);
        for (int j = 0; j < p; ++j) {
            beta_true[j] = truth[j].beta_star;
            xi_true[j]   = truth[j].xi;
        }
        K_true = 1 + *std::max_element(xi_true.begin(), xi_true.end());
        double s2 = 0.0;
        for (const auto& r : truth) s2 += r.u * r.u;
        tau_true = std::sqrt(s2 / p);
    }

    std::cerr << (truth_ok ? "Mode: SIMULATION (truth file provided)\n"
                           : "Mode: REAL DATA  (no ground truth)\n");
    std::cerr << "p=" << p << "  S=" << S << "  burn=" << S_burn
              << "  post-burnin=" << S_post << "\n\n";

    // =========================================================================
    // SECTION A — Scalar parameters
    // =========================================================================
    std::cout << "================================================================\n"
              << "SECTION A — Scalar parameter summaries\n"
              << "================================================================\n\n";

    auto print_scalar = [&](const std::string& name,
                            const std::vector<double>& chain,
                            bool has_tv, double tv = 0.0) -> Summary {
        std::vector<double> post(chain.begin() + S_burn, chain.end());
        const Summary s = summarise(post);
        const bool cov = has_tv && (s.lo95 <= tv && tv <= s.hi95);
        std::cout << std::left  << std::setw(8) << name
                  << "  mean=" << std::right << std::fixed << std::setprecision(5)
                  << std::setw(10) << s.mean
                  << "  sd=" << std::setw(8) << s.sd
                  << "  95%CI=[" << std::setw(10) << s.lo95
                  << ", " << std::setw(10) << s.hi95 << "]";
        if (has_tv) std::cout << "  truth=" << std::setw(10) << tv
                              << "  in_CI=" << (cov ? "YES" : "NO ");
        std::cout << "\n";
        return s;
    };

    print_scalar("gamma",  gamma_post, truth_ok, 0.0);
    print_scalar("psi",    psi_post,   false);
    print_scalar("tau",    tau_post,   truth_ok, tau_true);
    print_scalar("alpha",  alpha_post, false);

    // K
    {
        std::vector<double> Kc(S_post);
        for (int s = 0; s < S_post; ++s)
            Kc[s] = static_cast<double>(
                *std::max_element(xi_post[S_burn+s].begin(),
                                  xi_post[S_burn+s].end()) + 1);
        const Summary ks = summarise(Kc);
        std::cout << std::left << std::setw(8) << "K"
                  << "  mean=" << std::right << std::fixed << std::setprecision(5)
                  << std::setw(10) << ks.mean
                  << "  sd=" << std::setw(8) << ks.sd
                  << "  95%CI=[" << std::setw(10) << ks.lo95
                  << ", " << std::setw(10) << ks.hi95 << "]";
        if (truth_ok) std::cout << "  truth=" << std::setw(10) << K_true
                                << "  in_CI="
                                << (ks.lo95 <= K_true && K_true <= ks.hi95 ? "YES" : "NO ");
        std::cout << "\n";

        std::map<int,int> hist;
        for (double kd : Kc) ++hist[static_cast<int>(kd)];
        std::cout << "\n  Posterior distribution of K:\n";
        for (auto& [k, cnt] : hist) {
            const double pct = 100.0 * cnt / S_post;
            std::cout << "    K=" << std::setw(3) << k
                      << "  " << std::setw(6) << std::fixed << std::setprecision(2)
                      << pct << "%  |" << std::string(static_cast<int>(pct/2), '#') << "\n";
        }
    }

    // =========================================================================
    // SECTION B — Per-SNP causal effect beta
    // =========================================================================
    std::cout << "\n================================================================\n"
              << "SECTION B — Per-SNP causal effect beta\n"
              << "================================================================\n\n"
              << "  beta_hat[s][j] = beta_post[s][ xi_post[s][j] ]"
                 "  (label-switching free)\n\n";

    std::vector<double> bpm(p), blo(p), bhi(p);
    std::vector<bool>   bcov(p, false);

    for (int j = 0; j < p; ++j) {
        std::vector<double> draws(S_post);
        for (int s = 0; s < S_post; ++s)
            draws[s] = beta_post[S_burn+s][ xi_post[S_burn+s][j] ];
        const Summary sm = summarise(draws);
        bpm[j] = sm.mean; blo[j] = sm.lo95; bhi[j] = sm.hi95;
        if (truth_ok) bcov[j] = (sm.lo95 <= beta_true[j] && beta_true[j] <= sm.hi95);
    }

    if (truth_ok) {
        double rmse = 0.0; int ncov = 0;
        for (int j = 0; j < p; ++j) {
            rmse += (bpm[j]-beta_true[j])*(bpm[j]-beta_true[j]);
            if (bcov[j]) ++ncov;
        }
        std::cout << "  RMSE(beta_hat, beta_true) = "
                  << std::fixed << std::setprecision(6) << std::sqrt(rmse/p) << "\n"
                  << "  95% CI coverage           = " << std::setprecision(1)
                  << 100.0*ncov/p << "%  (" << ncov << "/" << p
                  << " SNPs)   [nominal: 95%]\n\n";
    }

    const int show = std::min(p, 15);
    std::cout << "  " << std::setw(5) << "SNP"
              << std::setw(12) << "post_mean"
              << std::setw(12) << "CI_lo95"
              << std::setw(12) << "CI_hi95";
    if (truth_ok) std::cout << std::setw(12) << "beta_true" << std::setw(9) << "in_CI";
    std::cout << "\n";
    for (int j = 0; j < show; ++j) {
        std::cout << "  " << std::setw(5) << (j+1)
                  << std::fixed << std::setprecision(5)
                  << std::setw(12) << bpm[j]
                  << std::setw(12) << blo[j]
                  << std::setw(12) << bhi[j];
        if (truth_ok) std::cout << std::setw(12) << beta_true[j]
                                << std::setw(9)  << (bcov[j] ? "YES" : "NO");
        std::cout << "\n";
    }
    if (show < p) std::cout << "  ... (" << p-show << " more SNPs in output file)\n";

    // =========================================================================
    // SECTION C — PSM and medoid partition
    // =========================================================================
    std::cout << "\n================================================================\n"
              << "SECTION C — Clustering  (PSM + medoid partition)\n"
              << "================================================================\n\n";

    std::cerr << "Computing PSM (" << p << "x" << p
              << ", " << S_post << " samples) ...\n";
    const auto psm = compute_psm(xi_post, S_burn, p);

    // PSM quality vs truth
    if (truth_ok) {
        double ps = 0.0, pd = 0.0; long ns = 0, nd = 0;
        for (int i = 0; i < p; ++i)
            for (int j = i+1; j < p; ++j) {
                const double v = psm[i*p+j];
                if (xi_true[i] == xi_true[j]) { ps += v; ++ns; }
                else                           { pd += v; ++nd; }
            }
        if (ns) ps /= ns;
        if (nd) pd /= nd;
        std::cout << "  PSM vs truth:\n"
                  << "    Mean PSM  same cluster    : " << std::fixed << std::setprecision(4)
                  << ps << "  [ideal 1.00]\n"
                  << "    Mean PSM  diff cluster    : " << pd << "  [ideal 0.00]\n"
                  << "    Discrimination gap        : " << ps-pd << "\n\n";
    }

    // PSM distribution summary
    {
        std::vector<double> od;
        od.reserve(p*(p-1)/2);
        for (int i = 0; i < p; ++i)
            for (int j = i+1; j < p; ++j) od.push_back(psm[i*p+j]);
        std::sort(od.begin(), od.end());
        const int np = static_cast<int>(od.size());
        int c90 = 0, c50 = 0, c10 = 0;
        for (double v : od) {
            if (v > 0.9) ++c90;
            if (v > 0.5) ++c50;
            if (v < 0.1) ++c10;
        }
        std::cout << "  PSM distribution (off-diagonal pairs):\n"
                  << "    PSM > 0.9 : " << std::fixed << std::setprecision(1)
                  << 100.0*c90/np << "%  (strongly co-clustered)\n"
                  << "    PSM > 0.5 : " << 100.0*c50/np << "%\n"
                  << "    PSM < 0.1 : " << 100.0*c10/np << "%  (strongly separated)\n\n";
    }

    // Medoid
    std::cerr << "Finding medoid partition ...\n";
    const auto xi_medoid = medoid_partition(xi_post, psm, S_burn, p);
    const int  K_med     = 1 + *std::max_element(xi_medoid.begin(), xi_medoid.end());

    std::cout << "  Medoid partition:  K_medoid = " << K_med;
    if (truth_ok) {
        std::cout << "  (K_true = " << K_true << ")\n"
                  << "  ARI(medoid, truth) = " << std::fixed << std::setprecision(4)
                  << adjusted_rand_index(xi_medoid, xi_true)
                  << "  [1.0 = perfect, ~0.0 = chance]\n";
    } else {
        std::cout << "\n";
    }

    // Cluster-level summary
    {
        std::map<int, std::vector<int>> meds;
        for (int j = 0; j < p; ++j) meds[xi_medoid[j]].push_back(j);
        std::cout << "\n  Medoid cluster summary:\n"
                  << "    " << std::setw(8) << "cluster"
                  << std::setw(7)  << "n_snps"
                  << std::setw(14) << "mean_beta_hat"
                  << std::setw(13) << "sd_beta_hat";
        if (truth_ok) std::cout << std::setw(15) << "mean_beta_true";
        std::cout << "\n";
        for (auto& [k, mems] : meds) {
            double mb = 0.0, mb2 = 0.0, mt = 0.0;
            for (int j : mems) { mb += bpm[j]; mb2 += bpm[j]*bpm[j]; if (truth_ok) mt += beta_true[j]; }
            const int sz = static_cast<int>(mems.size());
            mb /= sz; mb2 /= sz;
            const double sd = (sz > 1) ? std::sqrt(std::max(0.0, mb2 - mb*mb) * sz/(sz-1)) : 0.0;
            if (truth_ok) mt /= sz;
            std::cout << "    " << std::setw(8) << (k+1)
                      << std::setw(7) << sz
                      << std::fixed << std::setprecision(5)
                      << std::setw(14) << mb
                      << std::setw(13) << sd;
            if (truth_ok) std::cout << std::setw(15) << mt;
            std::cout << "\n";
        }
    }

    // =========================================================================
    // SECTION D — Binder loss trace
    // =========================================================================
    std::cout << "\n================================================================\n"
              << "SECTION D — Binder loss trace  (convergence check)\n"
              << "================================================================\n\n";

    std::vector<double> binder(S_post);
    double bmin = 1e300, bmax = 0.0;
    for (int s = 0; s < S_post; ++s) {
        const auto& xi = xi_post[S_burn+s];
        double d = 0.0;
        for (int i = 0; i < p; ++i)
            for (int j = i+1; j < p; ++j) {
                double diff = psm[i*p+j] - (xi[i]==xi[j] ? 1.0 : 0.0);
                d += diff * diff;
            }
        binder[s] = d;
        bmin = std::min(bmin, d);
        bmax = std::max(bmax, d);
    }

    // Running-mean bar chart
    std::cout << "  Running mean of Binder loss (should stabilise quickly after burn-in):\n";
    const int nr = std::min(S_post, 20);
    const int st = S_post / nr;
    double rsum = 0.0;
    for (int r = 0; r < nr; ++r) {
        const int end = (r+1)*st;
        for (int s = r*st; s < end; ++s) rsum += binder[s];
        const double rm = rsum / end;
        const double range = bmax - bmin;
        const int bar = static_cast<int>(40.0 * (1.0 - (range > 0 ? (rm-bmin)/range : 0.5)));
        std::cout << "    after iter " << std::setw(6) << (S_burn+end)
                  << "  mean=" << std::fixed << std::setprecision(2) << std::setw(10) << rm
                  << "  |" << std::string(bar,'=') << std::string(40-bar,' ') << "|\n";
    }

    // Trend: first-half vs second-half mean
    {
        double h1 = 0.0, h2 = 0.0;
        const int half = S_post / 2;
        for (int s = 0;    s < half;   ++s) h1 += binder[s];
        for (int s = half; s < S_post; ++s) h2 += binder[s];
        h1 /= half; h2 /= (S_post - half);
        const double pct = 100.0*(h2-h1) / (h1 > 0 ? h1 : 1.0);
        std::cout << "\n  First-half vs second-half mean: "
                  << std::fixed << std::setprecision(2) << h1 << " -> " << h2
                  << "  (" << std::showpos << std::setprecision(1) << pct << std::noshowpos << "%)\n"
                  << "  " << (std::abs(pct) < 10.0
                       ? "OK: no substantial trend — chain appears to have mixed."
                       : "WARNING: >10% trend — consider longer burn-in or more iterations.")
                  << "\n";
    }

    // =========================================================================
    // Write output files
    // =========================================================================
    std::cout << "\n================================================================\n"
              << "Writing output files ...\n"
              << "================================================================\n\n";

    // scalar CSV
    {
        std::ofstream f(opfx + "_scalar.csv");
        f << "param,mean,sd,ci_lo95,ci_hi95";
        if (truth_ok) f << ",truth,in_ci";
        f << "\n";
        auto wrow = [&](const std::string& nm,
                        const std::vector<double>& chain,
                        bool htv, double tv = 0.0) {
            std::vector<double> post(chain.begin()+S_burn, chain.end());
            const Summary s = summarise(post);
            f << nm << "," << s.mean << "," << s.sd << "," << s.lo95 << "," << s.hi95;
            if (truth_ok)
                f << "," << (htv ? std::to_string(tv) : "NA")
                  << "," << (htv ? ((s.lo95<=tv && tv<=s.hi95) ? "1" : "0") : "NA");
            f << "\n";
        };
        wrow("gamma", gamma_post, truth_ok, 0.0);
        wrow("psi",   psi_post,   false);
        wrow("tau",   tau_post,   truth_ok, tau_true);
        wrow("alpha", alpha_post, false);
        std::vector<double> Kc(S_post);
        for (int s = 0; s < S_post; ++s)
            Kc[s] = static_cast<double>(
                *std::max_element(xi_post[S_burn+s].begin(),
                                  xi_post[S_burn+s].end()) + 1);
        const Summary ks = summarise(Kc);
        f << "K," << ks.mean << "," << ks.sd << "," << ks.lo95 << "," << ks.hi95;
        if (truth_ok)
            f << "," << K_true << ","
              << (ks.lo95<=K_true && K_true<=ks.hi95 ? "1" : "0");
        f << "\n";
    }

    // beta_snp CSV
    {
        std::ofstream f(opfx + "_beta_snp.csv");
        f << "snp,beta_post_mean,ci_lo95,ci_hi95";
        if (truth_ok) f << ",xi_true,beta_true,in_ci";
        f << "\n";
        for (int j = 0; j < p; ++j) {
            f << (j+1) << "," << bpm[j] << "," << blo[j] << "," << bhi[j];
            if (truth_ok)
                f << "," << (xi_true[j]+1) << "," << beta_true[j]
                  << "," << (bcov[j] ? 1 : 0);
            f << "\n";
        }
    }

    // PSM (space-separated)
    {
        std::ofstream f(opfx + "_psm.csv");
        f << std::fixed << std::setprecision(6);
        for (int i = 0; i < p; ++i) {
            for (int j = 0; j < p; ++j) {
                f << psm[i*p+j];
                if (j+1 < p) f << " ";
            }
            f << "\n";
        }
    }

    // clustering CSV
    {
        std::ofstream f(opfx + "_clustering.csv");
        f << "snp,xi_medoid";
        if (truth_ok) f << ",xi_true,correct";
        f << "\n";
        std::map<int,int> bmap;
        if (truth_ok) {
            std::map<int, std::map<int,int>> cont;
            for (int j = 0; j < p; ++j) ++cont[xi_medoid[j]][xi_true[j]];
            for (auto& [km, inner] : cont) {
                int bt=-1, bn=-1;
                for (auto& [kt, n] : inner) if (n>bn) { bn=n; bt=kt; }
                bmap[km] = bt;
            }
        }
        for (int j = 0; j < p; ++j) {
            f << (j+1) << "," << (xi_medoid[j]+1);
            if (truth_ok)
                f << "," << (xi_true[j]+1) << ","
                  << (bmap[xi_medoid[j]] == xi_true[j] ? 1 : 0);
            f << "\n";
        }
    }

    // Binder loss trace
    {
        std::ofstream f(opfx + "_binder.csv");
        f << "iter,binder_loss\n";
        for (int s = 0; s < S_post; ++s)
            f << (S_burn+s+1) << "," << binder[s] << "\n";
    }

    std::cout << "  " << opfx << "_scalar.csv\n"
              << "  " << opfx << "_beta_snp.csv\n"
              << "  " << opfx << "_psm.csv      (" << p << "x" << p << " matrix)\n"
              << "  " << opfx << "_clustering.csv\n"
              << "  " << opfx << "_binder.csv\n\n";

    return 0;
}
