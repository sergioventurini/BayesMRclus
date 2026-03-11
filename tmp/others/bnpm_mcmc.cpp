// bnpm_mcmc.cpp
//
// MCMC sampler for Bayesian Nonparametric Model for Mendelian Randomization
// Models:
//   BNPM-MR  : DP mixture of causal effects, no pleiotropy (tau = 0)
//   BNPM-PMR : BNPM-MR + systematic pleiotropy overdispersion tau^2
//
// Algorithm components:
//   1. Exact Gibbs for latent instrument strengths gamma_j (per SNP)
//   2. NUTS (Hoffman & Gelman 2014) with dual averaging for (log psi, log tau)
//   3. Neal's Algorithm 8 for DP partition xi
//   4. Elliptical Slice Sampling (Murray et al. 2010) for cluster atoms beta_k*
//   5. Jain-Neal (2004) split-merge transdimensional moves
//   6. Escobar-West (1995) exact update for DP concentration alpha
//
// Compile: g++ -O2 -std=c++17 -o bnpm_mcmc bnpm_mcmc.cpp
// =============================================================================

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <vector>

// =============================================================================
// Section 1: Constants, global RNG, basic samplers
// =============================================================================

static constexpr double kLog2Pi   = 1.8378770664093454836;
static constexpr double kDeltaMax = 1000.0;   // NUTS divergence threshold

static std::mt19937_64 g_rng;

static inline double u01() {
    static std::uniform_real_distribution<double> d(0.0, 1.0);
    return d(g_rng);
}
static inline double stdnorm() {
    static std::normal_distribution<double> d(0.0, 1.0);
    return d(g_rng);
}
// Gamma(shape, scale=1/rate)
static inline double rgamma_dist(double shape, double scale) {
    std::gamma_distribution<double> d(shape, scale);
    return d(g_rng);
}
static inline double rbeta_dist(double a, double b) {
    double x = rgamma_dist(a, 1.0), y = rgamma_dist(b, 1.0);
    return x / (x + y);
}

// log-sum-exp
static double logsumexp(const std::vector<double>& v) {
    if (v.empty()) return -1e300;
    double m = *std::max_element(v.begin(), v.end());
    double s = 0.0;
    for (double x : v) s += std::exp(x - m);
    return m + std::log(s);
}

// Sample index from log-unnormalized probability vector
static int sample_lp(const std::vector<double>& lp) {
    double mx = *std::max_element(lp.begin(), lp.end());
    std::vector<double> p(lp.size());
    double sum = 0.0;
    for (int i = 0; i < (int)lp.size(); ++i) { p[i] = std::exp(lp[i] - mx); sum += p[i]; }
    double r = u01() * sum, c = 0.0;
    for (int i = 0; i < (int)lp.size(); ++i) {
        c += p[i];
        if (r <= c) return i;
    }
    return (int)lp.size() - 1;
}

// Dot product
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (int i = 0; i < (int)a.size(); ++i) s += a[i] * b[i];
    return s;
}

// =============================================================================
// Section 2: Data structures
// =============================================================================

struct SNP {
    double ghat;   // gamma-hat_j  (SNP-exposure effect estimate)
    double Ghat;   // Gamma-hat_j  (SNP-outcome effect estimate)
    double sg2;    // s^2_{gamma-hat_j} (squared SE for exposure)
    double sG2;    // s^2_{Gamma-hat_j} (squared SE for outcome)
};

struct HyperParams {
    double nu_psi   = 5.0;   // half-t df for psi
    double sig_psi  = 1.0;   // half-t scale for psi
    double nu_tau   = 5.0;   // half-t df for tau
    double sig_tau  = 1.0;   // half-t scale for tau
    double mu_beta  = 0.0;   // DP base measure mean
    double sig_beta = 1.0;   // DP base measure SD
    double a_alpha  = 1.0;   // Gamma prior shape for alpha
    double b_alpha  = 1.0;   // Gamma prior rate  for alpha
};

struct Config {
    int    S_wu         = 10000;  // warm-up iterations
    int    S            = 40000;  // sampling iterations
    int    T_SM         = 5;      // split-merge every T_SM sweeps
    int    m_aux        = 3;      // auxiliary clusters in Neal's Alg 8
    int    t_lnch       = 5;      // launch steps for split-merge
    double delta_nuts   = 0.65;   // NUTS target acceptance rate
    bool   use_tau      = true;   // false = BNPM-MR (tau fixed to 0)
    int    max_depth    = 10;     // NUTS maximum tree depth
    unsigned seed       = 42;
    std::string outfile = "bnpm_samples.csv";
};

// =============================================================================
// Section 3: Marginal log-likelihood (gamma_j integrated out analytically)
//
//  Model: gamma_j ~ N(0, psi^2)
//         gamma_hat_j = gamma_j + eps_j,  eps_j ~ N(0, sg2_j)
//         Gamma_hat_j = b_j*gamma_j + eta_j, eta_j ~ N(0, sG2_j + tau^2)
//
//  Marginalising over gamma_j:
//    [gamma_hat_j, Gamma_hat_j]^T | beta_j, psi, tau ~ N_2(0, Sigma_j)
//    where Sigma_j = [[a_j, c_j], [c_j, d_j]]
//          a_j = psi^2 + sg2_j
//          c_j = b * psi^2
//          d_j = b^2 * psi^2 + sG2_j + tau^2
//          D_j = a_j*d_j - c_j^2   (determinant)
//    log_lik_j = -log(2*pi) - 0.5*log(D_j)
//                - 0.5*(d_j*u^2 - 2*c_j*u*v + a_j*v^2)/D_j
//    u = gamma_hat_j, v = Gamma_hat_j
// =============================================================================

static double loglik1(const SNP& s, double b, double psi, double tau) {
    double u   = s.ghat, v = s.Ghat;
    double psi2 = psi * psi, tau2 = tau * tau;
    double a = psi2 + s.sg2;
    double c = b * psi2;
    double d = b * b * psi2 + s.sG2 + tau2;
    double D = a * d - c * c;
    if (D <= 0.0) return -1e300;
    double Q = (d * u * u - 2.0 * c * u * v + a * v * v) / D;
    return -kLog2Pi - 0.5 * std::log(D) - 0.5 * Q;
}

// Sum of log-likelihoods over a set of SNP indices
static double loglik_set(const std::vector<SNP>& data,
                          const std::vector<int>& idx,
                          double b, double psi, double tau) {
    double ll = 0.0;
    for (int j : idx) ll += loglik1(data[j], b, psi, tau);
    return ll;
}

// =============================================================================
// Section 4: Gradients of marginal log-likelihood
//   w.r.t. psi^2 and tau^2 (used by NUTS via chain rule)
// =============================================================================

static void grad_loglik1(const SNP& s, double b, double psi, double tau,
                          double& dl_dpsi2, double& dl_dtau2) {
    double u   = s.ghat, v = s.Ghat;
    double psi2 = psi * psi, tau2 = tau * tau;
    double a = psi2 + s.sg2;
    double c = b * psi2;
    double d = b * b * psi2 + s.sG2 + tau2;
    double D = a * d - c * c;
    if (D <= 0.0) { dl_dpsi2 = 0.0; dl_dtau2 = 0.0; return; }
    double Q  = (d * u * u - 2.0 * c * u * v + a * v * v) / D;

    // dD/d(psi^2) = d + a*b^2 - 2*c*b = sG2 + tau^2 + b^2*sg2
    double dD_dpsi2 = s.sG2 + tau2 + b * b * s.sg2;
    // dN/d(psi^2) = (b*u - v)^2   (only d,c,a depend on psi^2 in N)
    double dN_dpsi2 = (b * u - v) * (b * u - v);
    double dQ_dpsi2 = dN_dpsi2 / D - Q * dD_dpsi2 / D;
    dl_dpsi2 = -0.5 * dD_dpsi2 / D - 0.5 * dQ_dpsi2;

    // dD/d(tau^2) = a_j
    // dN/d(tau^2) = u^2  (only d_j = b^2*psi^2 + sG2 + tau^2 changes)
    double dN_dtau2 = u * u;
    double dQ_dtau2 = dN_dtau2 / D - Q * a / D;
    dl_dtau2 = -0.5 * a / D - 0.5 * dQ_dtau2;
}

// =============================================================================
// Section 5: Potential energy U and gradient for NUTS
//
//  theta[0] = log(psi);  theta[1] = log(tau)  [if use_tau]
//
//  U(theta) = -sum_j ll_j(theta) - log p(theta)
//
//  Half-t prior: p(psi) propto (1 + psi^2/(nu*sig^2))^{-(nu+1)/2}  for psi>0
//  Reparameterised to log_psi = theta[0]:
//    log p(theta[0]) = theta[0] - (nu+1)/2 * log(1 + e^{2*theta[0]}/(nu*sig^2))
//    d/d(theta[0]) = 1 - (nu+1)*psi^2 / (nu*sig^2 + psi^2)
// =============================================================================

static double potential(const std::vector<double>& theta,
                         const std::vector<SNP>& data,
                         const std::vector<int>& xi,
                         const std::vector<double>& atoms,
                         const HyperParams& hp,
                         bool use_tau) {
    double psi = std::exp(theta[0]);
    double tau = (use_tau && (int)theta.size() > 1) ? std::exp(theta[1]) : 0.0;

    double U = 0.0;
    for (int j = 0; j < (int)data.size(); ++j)
        U -= loglik1(data[j], atoms[xi[j]], psi, tau);

    // Half-t prior on psi (reparameterised to log scale)
    double psi2 = psi * psi;
    double nu_p = hp.nu_psi, s2_p = hp.sig_psi * hp.sig_psi;
    U += 0.5 * (nu_p + 1.0) * std::log(1.0 + psi2 / (nu_p * s2_p)) - theta[0];

    if (use_tau && (int)theta.size() > 1) {
        double tau2 = tau * tau;
        double nu_t = hp.nu_tau, s2_t = hp.sig_tau * hp.sig_tau;
        U += 0.5 * (nu_t + 1.0) * std::log(1.0 + tau2 / (nu_t * s2_t)) - theta[1];
    }
    return U;
}

static std::vector<double> grad_potential(const std::vector<double>& theta,
                                            const std::vector<SNP>& data,
                                            const std::vector<int>& xi,
                                            const std::vector<double>& atoms,
                                            const HyperParams& hp,
                                            bool use_tau) {
    int dim = (int)theta.size();
    std::vector<double> grad(dim, 0.0);
    double psi  = std::exp(theta[0]);
    double tau  = (use_tau && dim > 1) ? std::exp(theta[1]) : 0.0;
    double psi2 = psi * psi;

    double dU_dpsi2 = 0.0, dU_dtau2 = 0.0;
    for (int j = 0; j < (int)data.size(); ++j) {
        double dl_dpsi2, dl_dtau2;
        grad_loglik1(data[j], atoms[xi[j]], psi, tau, dl_dpsi2, dl_dtau2);
        dU_dpsi2 -= dl_dpsi2;
        dU_dtau2 -= dl_dtau2;
    }

    // Chain rule: dU/d(log_psi) = dU/d(psi^2) * 2*psi^2
    grad[0] = dU_dpsi2 * 2.0 * psi2;
    // Add prior gradient for log_psi
    double nu_p = hp.nu_psi, s2_p = hp.sig_psi * hp.sig_psi;
    grad[0] += (nu_p + 1.0) * psi2 / (nu_p * s2_p + psi2) - 1.0;

    if (use_tau && dim > 1) {
        double tau2 = tau * tau;
        grad[1] = dU_dtau2 * 2.0 * tau2;
        double nu_t = hp.nu_tau, s2_t = hp.sig_tau * hp.sig_tau;
        grad[1] += (nu_t + 1.0) * tau2 / (nu_t * s2_t + tau2) - 1.0;
    }
    return grad;
}

// =============================================================================
// Section 6: NUTS (No-U-Turn Sampler) with dual averaging
//   Implements Algorithm 3 + 6 of Hoffman & Gelman (2014)
// =============================================================================

static void leapfrog(std::vector<double>& theta,
                      std::vector<double>& r,
                      double eps,
                      const std::vector<SNP>& data,
                      const std::vector<int>& xi,
                      const std::vector<double>& atoms,
                      const HyperParams& hp,
                      bool use_tau) {
    auto g = grad_potential(theta, data, xi, atoms, hp, use_tau);
    for (int i = 0; i < (int)r.size(); ++i) r[i]     -= 0.5 * eps * g[i];
    for (int i = 0; i < (int)theta.size(); ++i) theta[i] += eps * r[i];
    g = grad_potential(theta, data, xi, atoms, hp, use_tau);
    for (int i = 0; i < (int)r.size(); ++i) r[i]     -= 0.5 * eps * g[i];
}

static double hamiltonian(const std::vector<double>& theta,
                            const std::vector<double>& r,
                            const std::vector<SNP>& data,
                            const std::vector<int>& xi,
                            const std::vector<double>& atoms,
                            const HyperParams& hp,
                            bool use_tau) {
    return potential(theta, data, xi, atoms, hp, use_tau) + 0.5 * dot(r, r);
}

// Result of one build_tree call
struct BTRes {
    std::vector<double> theta_m, theta_p;  // minus/plus trajectory endpoints
    std::vector<double> r_m,     r_p;
    std::vector<double> theta_prop;        // proposed sample
    int    n;           // number of valid (in-slice) states
    bool   s;           // stop criterion: no U-turn yet
    double sum_alpha;   // cumulative MH acceptance prob (for adaptation)
    int    n_alpha;     // count of states in sum_alpha
};

// Recursive NUTS tree builder (Hoffman & Gelman 2014, Algorithm 6)
static BTRes build_tree(std::vector<double> theta,
                         std::vector<double> r,
                         double log_u,
                         int v,        // direction: -1 (backward) or +1 (forward)
                         int j,        // tree depth
                         double eps,
                         double H0,
                         const std::vector<SNP>& data,
                         const std::vector<int>& xi,
                         const std::vector<double>& atoms,
                         const HyperParams& hp,
                         bool use_tau) {
    if (j == 0) {
        // Base case: single leapfrog step in direction v
        leapfrog(theta, r, v * eps, data, xi, atoms, hp, use_tau);
        double H_new = hamiltonian(theta, r, data, xi, atoms, hp, use_tau);
        double dH    = H0 - H_new;
        BTRes res;
        res.theta_m    = res.theta_p    = theta;
        res.r_m        = res.r_p        = r;
        res.theta_prop = theta;
        res.n          = (log_u <= -H_new) ? 1 : 0;
        res.s          = (log_u < kDeltaMax - H_new);
        res.sum_alpha  = std::min(1.0, std::exp(dH));
        res.n_alpha    = 1;
        return res;
    }

    // Recursive case: build first half-subtree
    BTRes res1 = build_tree(theta, r, log_u, v, j - 1, eps, H0,
                             data, xi, atoms, hp, use_tau);
    if (!res1.s) return res1;  // early stopping

    // Build second half-subtree from appropriate endpoint
    BTRes res2;
    if (v == -1)
        res2 = build_tree(res1.theta_m, res1.r_m, log_u, v, j - 1, eps, H0,
                           data, xi, atoms, hp, use_tau);
    else
        res2 = build_tree(res1.theta_p, res1.r_p, log_u, v, j - 1, eps, H0,
                           data, xi, atoms, hp, use_tau);

    // Progressive Metropolis: accept res2's proposal with prob n2/(n1+n2)
    int n_tot = res1.n + res2.n;
    BTRes res;
    if (n_tot > 0 && u01() < (double)res2.n / (double)n_tot)
        res.theta_prop = res2.theta_prop;
    else
        res.theta_prop = res1.theta_prop;

    // Set trajectory endpoints based on direction
    if (v == -1) {
        res.theta_m = res2.theta_m; res.r_m = res2.r_m;
        res.theta_p = res1.theta_p; res.r_p = res1.r_p;
    } else {
        res.theta_m = res1.theta_m; res.r_m = res1.r_m;
        res.theta_p = res2.theta_p; res.r_p = res2.r_p;
    }

    // U-turn check: (theta+ - theta-) . r- >= 0  AND  (theta+ - theta-) . r+ >= 0
    int dim = (int)theta.size();
    std::vector<double> dth(dim);
    for (int i = 0; i < dim; ++i) dth[i] = res.theta_p[i] - res.theta_m[i];
    bool uturn = (dot(dth, res.r_m) < 0.0) || (dot(dth, res.r_p) < 0.0);

    res.n       = n_tot;
    res.s       = res1.s && res2.s && !uturn;
    res.sum_alpha = res1.sum_alpha + res2.sum_alpha;
    res.n_alpha   = res1.n_alpha   + res2.n_alpha;
    return res;
}

// Full NUTS transition; returns mean acceptance probability for dual averaging
static double nuts_step(std::vector<double>& theta,
                         double eps,
                         const std::vector<SNP>& data,
                         const std::vector<int>& xi,
                         const std::vector<double>& atoms,
                         const HyperParams& hp,
                         bool use_tau,
                         int max_depth) {
    int dim = (int)theta.size();
    std::vector<double> r(dim);
    for (int i = 0; i < dim; ++i) r[i] = stdnorm();

    double H0    = hamiltonian(theta, r, data, xi, atoms, hp, use_tau);
    double log_u = std::log(u01()) - H0;

    std::vector<double> theta_m = theta, theta_p = theta;
    std::vector<double> r_m     = r,     r_p     = r;
    std::vector<double> theta_cur = theta;
    int    n          = 1;
    bool   s          = true;
    double sum_alpha  = 0.0;
    int    n_alpha    = 0;

    for (int j = 0; j < max_depth && s; ++j) {
        int v = (u01() < 0.5) ? -1 : 1;
        BTRes res;
        if (v == -1)
            res = build_tree(theta_m, r_m, log_u, v, j, eps, H0,
                              data, xi, atoms, hp, use_tau);
        else
            res = build_tree(theta_p, r_p, log_u, v, j, eps, H0,
                              data, xi, atoms, hp, use_tau);

        sum_alpha += res.sum_alpha;
        n_alpha   += res.n_alpha;

        if (res.s && res.n > 0) {
            double accept = std::min(1.0, (double)res.n / (double)n);
            if (u01() < accept) theta_cur = res.theta_prop;
        }

        if (v == -1) { theta_m = res.theta_m; r_m = res.r_m; }
        else         { theta_p = res.theta_p; r_p = res.r_p; }
        n += res.n;

        // Global U-turn check at each level
        std::vector<double> dth(dim);
        for (int i = 0; i < dim; ++i) dth[i] = theta_p[i] - theta_m[i];
        s = res.s && (dot(dth, r_m) >= 0.0) && (dot(dth, r_p) >= 0.0);
    }
    theta = theta_cur;
    return (n_alpha > 0) ? sum_alpha / (double)n_alpha : 0.0;
}

// Nesterov dual averaging for step-size adaptation (Hoffman & Gelman 2014, Alg 5)
static void dual_avg_update(double& log_eps,
                              double& log_eps_bar,
                              double& H_bar,
                              int m,          // current iteration (1-indexed)
                              double mu,      // log(10 * eps0)
                              double alpha_bar,
                              double delta,
                              double gamma_da = 0.05,
                              double t0       = 10.0,
                              double kappa    = 0.75) {
    double md = (double)m;
    H_bar     = (1.0 - 1.0 / (md + t0)) * H_bar + (delta - alpha_bar) / (md + t0);
    log_eps   = mu - std::sqrt(md) / gamma_da * H_bar;
    double mk = std::pow(md, -kappa);
    log_eps_bar = mk * log_eps + (1.0 - mk) * log_eps_bar;
}

// =============================================================================
// Section 7: Exact Gibbs update for latent gamma_j
//
//  Model conditional on beta_j, psi, tau, gamma_hat_j, Gamma_hat_j:
//    p(gamma_j | ...) = N(mu_j, v_j^2)
//    1/v_j^2 = 1/psi^2 + 1/sg2_j + beta_j^2 / (sG2_j + tau^2)
//    mu_j    = v_j^2 * (gamma_hat_j/sg2_j + beta_j*Gamma_hat_j/(sG2_j+tau^2))
// =============================================================================

static void update_gamma(std::vector<double>& gamma,
                          const std::vector<SNP>& data,
                          const std::vector<int>& xi,
                          const std::vector<double>& atoms,
                          double psi, double tau) {
    double psi2 = psi * psi, tau2 = tau * tau;
    for (int j = 0; j < (int)data.size(); ++j) {
        double b      = atoms[xi[j]];
        double sG2eff = data[j].sG2 + tau2;
        double prec   = 1.0 / psi2 + 1.0 / data[j].sg2 + b * b / sG2eff;
        double mu     = (data[j].ghat / data[j].sg2 + b * data[j].Ghat / sG2eff) / prec;
        gamma[j]      = mu + stdnorm() / std::sqrt(prec);
    }
}

// =============================================================================
// Section 8: Cluster management utilities
// =============================================================================

static std::vector<std::vector<int>> cluster_members(const std::vector<int>& xi, int K) {
    std::vector<std::vector<int>> mem(K);
    for (int j = 0; j < (int)xi.size(); ++j)
        if (xi[j] >= 0 && xi[j] < K) mem[xi[j]].push_back(j);
    return mem;
}

// Remove empty clusters and relabel 0..K-1 (compact representation)
static void compact(std::vector<int>& xi, std::vector<double>& atoms, int& K) {
    std::vector<int> new_label(K, -1);
    int new_K = 0;
    std::vector<double> new_atoms;
    for (int k = 0; k < K; ++k) {
        bool occ = false;
        for (int j : xi) if (j == k) { occ = true; break; }
        if (occ) { new_label[k] = new_K++; new_atoms.push_back(atoms[k]); }
    }
    for (int& j : xi) { assert(new_label[j] >= 0); j = new_label[j]; }
    atoms = new_atoms;
    K     = new_K;
}

// =============================================================================
// Section 9: Neal's Algorithm 8 for DP partition update
//   Uses m_aux auxiliary tables drawn fresh from the base measure G_0
// =============================================================================

static void update_partition_neal8(std::vector<int>& xi,
                                     std::vector<double>& atoms,
                                     int& K,
                                     const std::vector<SNP>& data,
                                     double psi, double tau,
                                     double alpha,
                                     const HyperParams& hp,
                                     int m_aux) {
    int p = (int)data.size();

    for (int j = 0; j < p; ++j) {
        int k_old = xi[j];
        xi[j] = -1;  // temporarily remove j

        // Cluster counts with j removed
        std::vector<int> cnt(K, 0);
        for (int jj = 0; jj < p; ++jj)
            if (xi[jj] >= 0 && xi[jj] < K) cnt[xi[jj]]++;

        bool is_sing = (cnt[k_old] == 0);

        // Build candidate list: {atom, log_prob, existing cluster index or -1}
        struct Cand { double atom; double lp; int k; };
        std::vector<Cand> cands;

        // Existing non-empty clusters: weight = n_k * lik(j | beta_k*)
        for (int k = 0; k < K; ++k)
            if (cnt[k] > 0)
                cands.push_back({atoms[k],
                                  std::log((double)cnt[k]) + loglik1(data[j], atoms[k], psi, tau),
                                  k});

        // Auxiliary clusters: weight = (alpha/m) * lik(j | b_aux)
        if (is_sing) {
            // First auxiliary slot reuses the orphaned atom (standard in Alg 8)
            cands.push_back({atoms[k_old],
                              std::log(alpha / (double)m_aux) + loglik1(data[j], atoms[k_old], psi, tau),
                              -1});
            for (int t = 0; t < m_aux - 1; ++t) {
                double b = hp.mu_beta + hp.sig_beta * stdnorm();
                cands.push_back({b,
                                  std::log(alpha / (double)m_aux) + loglik1(data[j], b, psi, tau),
                                  -1});
            }
        } else {
            for (int t = 0; t < m_aux; ++t) {
                double b = hp.mu_beta + hp.sig_beta * stdnorm();
                cands.push_back({b,
                                  std::log(alpha / (double)m_aux) + loglik1(data[j], b, psi, tau),
                                  -1});
            }
        }

        // Sample a candidate
        std::vector<double> lps;
        for (auto& c : cands) lps.push_back(c.lp);
        int idx = sample_lp(lps);

        if (cands[idx].k >= 0) {
            // Assign to existing cluster
            xi[j] = cands[idx].k;
            // k_old may now be empty; compact() will handle it at end of sweep
        } else {
            // New cluster
            if (is_sing) {
                // Reuse k_old slot (avoids fragmenting indices)
                atoms[k_old] = cands[idx].atom;
                xi[j]        = k_old;
            } else {
                atoms.push_back(cands[idx].atom);
                xi[j] = K++;
            }
        }
    }
    compact(xi, atoms, K);
}

// =============================================================================
// Section 10: Elliptical Slice Sampling (Murray et al. 2010) for cluster atoms
//   Prior:        beta_k* ~ N(mu_beta, sig_beta^2)
//   Log-lik:      ell(b) = sum_{j in cluster k} loglik1(data[j], b, psi, tau)
//   ESS exploits the Gaussian prior by tracing an ellipse in parameter space
// =============================================================================

// ESS for a single atom; returns updated value
static double ess_one(double beta_cur,
                       const std::vector<int>& members,
                       const std::vector<SNP>& data,
                       double psi, double tau,
                       double mu_beta, double sig_beta) {
    // Current log-likelihood
    double ll_cur = 0.0;
    for (int j : members) ll_cur += loglik1(data[j], beta_cur, psi, tau);

    double nu    = mu_beta + sig_beta * stdnorm();  // prior sample (ellipse partner)
    double log_y = ll_cur + std::log(u01());         // slice threshold

    // Initial angle bracket
    double theta_ = (u01() * 2.0 - 1.0) * M_PI;    // Uniform(-pi, pi)
    double t_min  = theta_ - 2.0 * M_PI;
    double t_max  = theta_;

    double f0 = beta_cur - mu_beta;   // centered current
    double nc  = nu       - mu_beta;  // centered prior draw

    for (int iter = 0; iter < 200; ++iter) {
        // Proposal on the ellipse: beta' = f0*cos(theta) + nc*sin(theta) + mu_beta
        double bp  = f0 * std::cos(theta_) + nc * std::sin(theta_) + mu_beta;
        double ll  = 0.0;
        for (int j : members) ll += loglik1(data[j], bp, psi, tau);
        if (ll >= log_y) return bp;           // accepted
        // Shrink bracket
        if (theta_ < 0.0) t_min = theta_;
        else               t_max = theta_;
        theta_ = t_min + u01() * (t_max - t_min);
    }
    return beta_cur;  // fallback (extremely rare)
}

// ESS update for all cluster atoms
static void update_atoms_ess(std::vector<double>& atoms,
                               const std::vector<int>& xi,
                               int K,
                               const std::vector<SNP>& data,
                               double psi, double tau,
                               const HyperParams& hp) {
    auto mem = cluster_members(xi, K);
    for (int k = 0; k < K; ++k)
        if (!mem[k].empty())
            atoms[k] = ess_one(atoms[k], mem[k], data, psi, tau,
                                hp.mu_beta, hp.sig_beta);
}

// =============================================================================
// Section 11: Jain-Neal (2004) split-merge transdimensional moves
//
//  Select two random SNPs i and j:
//    If xi[i] == xi[j]:  propose SPLIT  (one cluster -> two)
//    Else:               propose MERGE  (two clusters -> one)
//
//  Uses a launch sequence of t_lnch restricted Gibbs sweeps + ESS atom updates
//  to improve proposal quality and compute the asymmetric proposal probabilities.
// =============================================================================

// Compute log probability of the current restricted allocation S
// given atoms beta_a, beta_b and final cluster counts.
// i_fixed is always in S=0; j_fixed is always in S=1.
static double log_q_allocation(const std::vector<SNP>& data,
                                 const std::vector<int>& local_idx,
                                 const std::vector<int>& S,   // 0 or 1
                                 int i_fixed, int j_fixed,
                                 double beta_a, double beta_b,
                                 double psi, double tau) {
    int n_a = 0, n_b = 0;
    for (int l : local_idx) { if (S[l] == 0) n_a++; else n_b++; }

    double lq = 0.0;
    for (int l : local_idx) {
        if (l == i_fixed || l == j_fixed) continue;
        int na_m = n_a - (S[l] == 0 ? 1 : 0);
        int nb_m = n_b - (S[l] == 1 ? 1 : 0);
        double lp_a = (na_m > 0 ? std::log((double)na_m) : -1e300)
                      + loglik1(data[l], beta_a, psi, tau);
        double lp_b = (nb_m > 0 ? std::log((double)nb_m) : -1e300)
                      + loglik1(data[l], beta_b, psi, tau);
        double lse  = logsumexp({lp_a, lp_b});
        lq += (S[l] == 0) ? (lp_a - lse) : (lp_b - lse);
    }
    return lq;
}

// One restricted Gibbs sweep (updates S in place; keeps i_fixed in 0, j_fixed in 1)
static void restricted_gibbs_sweep(const std::vector<SNP>& data,
                                     const std::vector<int>& local_idx,
                                     std::vector<int>& S,
                                     int i_fixed, int j_fixed,
                                     double beta_a, double beta_b,
                                     double psi, double tau,
                                     int& n_a, int& n_b) {
    for (int l : local_idx) {
        if (l == i_fixed || l == j_fixed) continue;
        int na_m = n_a - (S[l] == 0 ? 1 : 0);
        int nb_m = n_b - (S[l] == 1 ? 1 : 0);
        double lp_a = (na_m > 0 ? std::log((double)na_m) : -1e300)
                      + loglik1(data[l], beta_a, psi, tau);
        double lp_b = (nb_m > 0 ? std::log((double)nb_m) : -1e300)
                      + loglik1(data[l], beta_b, psi, tau);
        int old_s = S[l];
        int new_s = sample_lp({lp_a, lp_b});
        if (new_s != old_s) {
            if (new_s == 0) { n_a++; n_b--; } else { n_a--; n_b++; }
            S[l] = new_s;
        }
    }
}

// Run launch sequence: t_lnch restricted Gibbs sweeps interleaved with ESS atom updates
static void run_launch(const std::vector<SNP>& data,
                        const std::vector<int>& local_idx,
                        std::vector<int>& S,
                        int i_fixed, int j_fixed,
                        double& beta_a, double& beta_b,
                        double psi, double tau,
                        int t_lnch,
                        double mu_beta, double sig_beta,
                        int& n_a, int& n_b) {
    for (int t = 0; t < t_lnch; ++t) {
        restricted_gibbs_sweep(data, local_idx, S, i_fixed, j_fixed,
                                beta_a, beta_b, psi, tau, n_a, n_b);
        // ESS atom updates after each Gibbs sweep
        std::vector<int> mem_a, mem_b;
        for (int l : local_idx) {
            if      (S[l] == 0) mem_a.push_back(l);
            else if (S[l] == 1) mem_b.push_back(l);
        }
        if (!mem_a.empty()) beta_a = ess_one(beta_a, mem_a, data, psi, tau, mu_beta, sig_beta);
        if (!mem_b.empty()) beta_b = ess_one(beta_b, mem_b, data, psi, tau, mu_beta, sig_beta);
    }
}

// Main split-merge move
static void split_merge_move(std::vector<int>& xi,
                               std::vector<double>& atoms,
                               int& K,
                               const std::vector<SNP>& data,
                               double psi, double tau,
                               double alpha,
                               const HyperParams& hp,
                               int t_lnch) {
    int p = (int)data.size();
    if (p < 2) return;

    // Draw two distinct SNP indices
    int idx_i = (int)(u01() * p) % p;
    int idx_j;
    do { idx_j = (int)(u01() * p) % p; } while (idx_j == idx_i);

    int ki = xi[idx_i], kj = xi[idx_j];

    auto log_prior_beta = [&](double b) -> double {
        double z = (b - hp.mu_beta) / hp.sig_beta;
        return -0.5 * z * z - std::log(hp.sig_beta) - 0.5 * std::log(2.0 * M_PI);
    };

    if (ki == kj) {
        // ================================================================
        // SPLIT MOVE: propose splitting cluster ki into two
        // ================================================================
        int kc = ki;
        auto mem = cluster_members(xi, K);
        auto& cl = mem[kc];
        int n_c  = (int)cl.size();
        if (n_c < 2) return;

        // Draw new atoms from the prior
        double beta_a = hp.mu_beta + hp.sig_beta * stdnorm();
        double beta_b = hp.mu_beta + hp.sig_beta * stdnorm();

        // Random initial allocation (idx_i in S_a=0, idx_j in S_b=1)
        std::vector<int> S(p, -1);
        int n_a = 1, n_b = 1;
        S[idx_i] = 0; S[idx_j] = 1;
        for (int l : cl) {
            if (l == idx_i || l == idx_j) continue;
            if (u01() < 0.5) { S[l] = 0; n_a++; } else { S[l] = 1; n_b++; }
        }

        // Launch sequence
        run_launch(data, cl, S, idx_i, idx_j, beta_a, beta_b,
                   psi, tau, t_lnch, hp.mu_beta, hp.sig_beta, n_a, n_b);

        // Proposal probability: product of final conditional probabilities
        double lq_split = log_q_allocation(data, cl, S, idx_i, idx_j,
                                            beta_a, beta_b, psi, tau);

        // Gather final sub-cluster members
        std::vector<int> mem_a, mem_b;
        for (int l : cl) {
            if      (S[l] == 0) mem_a.push_back(l);
            else if (S[l] == 1) mem_b.push_back(l);
        }
        int ni = n_a, nj = n_b;

        // Log-likelihoods under proposed and current state
        double ll_a = loglik_set(data, mem_a, beta_a,     psi, tau);
        double ll_b = loglik_set(data, mem_b, beta_b,     psi, tau);
        double ll_c = loglik_set(data, cl,    atoms[kc],  psi, tau);

        // Log acceptance ratio (Jain-Neal 2004, Eq. 13)
        double log_A =
              std::lgamma((double)ni) + std::lgamma((double)nj) - std::lgamma((double)n_c)
            + std::log(alpha)
            + ll_a + ll_b - ll_c
            + log_prior_beta(beta_a) + log_prior_beta(beta_b) - log_prior_beta(atoms[kc])
            - lq_split;

        if (std::log(u01()) < log_A) {
            // Accept: reuse kc for sub-cluster a; add new cluster for b
            atoms[kc] = beta_a;
            atoms.push_back(beta_b);
            int k_new = K++;
            for (int l : cl) xi[l] = (S[l] == 0) ? kc : k_new;
        }

    } else {
        // ================================================================
        // MERGE MOVE: propose merging clusters ki and kj into one
        // ================================================================
        int kc = ki, kd = kj;
        auto mem = cluster_members(xi, K);
        auto& cl_c = mem[kc];
        auto& cl_d = mem[kd];
        int n_c = (int)cl_c.size();
        int n_d = (int)cl_d.size();
        int n_m = n_c + n_d;

        // Combined list
        std::vector<int> cl_cd;
        cl_cd.insert(cl_cd.end(), cl_c.begin(), cl_c.end());
        cl_cd.insert(cl_cd.end(), cl_d.begin(), cl_d.end());

        // Propose merge atom from prior
        double beta_new = hp.mu_beta + hp.sig_beta * stdnorm();

        // Compute reverse proposal probability (what would split look like?)
        // Start from original allocation and run launch to get q_split_reverse
        std::vector<int> S_rev(p, -1);
        for (int l : cl_c) S_rev[l] = 0;
        for (int l : cl_d) S_rev[l] = 1;
        int na_r = n_c, nb_r = n_d;
        double beta_a_rev = atoms[kc], beta_b_rev = atoms[kd];

        run_launch(data, cl_cd, S_rev, idx_i, idx_j, beta_a_rev, beta_b_rev,
                   psi, tau, t_lnch, hp.mu_beta, hp.sig_beta, na_r, nb_r);

        // Evaluate reverse probability at original allocation
        std::vector<int> S_orig(p, -1);
        for (int l : cl_c) S_orig[l] = 0;
        for (int l : cl_d) S_orig[l] = 1;
        double lq_rev = log_q_allocation(data, cl_cd, S_orig, idx_i, idx_j,
                                          atoms[kc], atoms[kd], psi, tau);

        // Log-likelihoods
        double ll_c   = loglik_set(data, cl_c,  atoms[kc],  psi, tau);
        double ll_d   = loglik_set(data, cl_d,  atoms[kd],  psi, tau);
        double ll_new = loglik_set(data, cl_cd, beta_new,   psi, tau);

        // Log acceptance ratio for merge (reverse of split)
        double log_A =
              std::lgamma((double)n_m) - std::lgamma((double)n_c) - std::lgamma((double)n_d)
            - std::log(alpha)
            + ll_new - ll_c - ll_d
            + log_prior_beta(beta_new) - log_prior_beta(atoms[kc]) - log_prior_beta(atoms[kd])
            + lq_rev;

        if (std::log(u01()) < log_A) {
            // Accept: merge kd into kc
            atoms[kc] = beta_new;
            for (int l : cl_d) xi[l] = kc;
            compact(xi, atoms, K);
        }
    }
}

// =============================================================================
// Section 12: Escobar-West (1995) exact update for DP concentration alpha
//   Auxiliary variable trick with exact two-component mixture
// =============================================================================

static void update_alpha(double& alpha, int K, int p, const HyperParams& hp) {
    double a   = hp.a_alpha, b = hp.b_alpha;
    double eta = rbeta_dist(alpha + 1.0, (double)p);
    double b_new = b - std::log(eta);
    // Mixture weight (see Escobar & West 1995, eq. A.2)
    double pi1 = (a + (double)K - 1.0) /
                 (a + (double)K - 1.0 + (double)p * b_new);
    if (u01() < pi1)
        alpha = rgamma_dist(a + (double)K,       1.0 / b_new);
    else
        alpha = rgamma_dist(a + (double)K - 1.0, 1.0 / b_new);
    alpha = std::max(alpha, 1e-6);
}

// =============================================================================
// Section 13: Synthetic data generator
// =============================================================================

static std::vector<SNP> generate_data(int p,
                                        double psi_true,
                                        double tau_true,
                                        const std::vector<double>& betas,
                                        const std::vector<double>& pi_true,
                                        unsigned seed) {
    std::mt19937_64 rng_loc(seed);
    std::normal_distribution<double> norm(0.0, 1.0);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::uniform_real_distribution<double> unif_se_exp(0.03, 0.07);
    std::uniform_real_distribution<double> unif_sG_exp(0.05, 0.09);

    std::vector<SNP> data(p);
    for (int j = 0; j < p; ++j) {
        // Component assignment via multinomial
        double r   = unif(rng_loc);
        double cum = 0.0;
        int k = (int)pi_true.size() - 1;
        for (int c = 0; c < (int)pi_true.size(); ++c) {
            cum += pi_true[c];
            if (r <= cum) { k = c; break; }
        }
        double b   = betas[k];
        double sg  = unif_se_exp(rng_loc);   // SE for exposure
        double sG  = unif_sG_exp(rng_loc);   // SE for outcome
        data[j].sg2 = sg * sg;
        data[j].sG2 = sG * sG;
        double gj   = psi_true * norm(rng_loc);  // true gamma_j
        data[j].ghat = gj + sg * norm(rng_loc);
        data[j].Ghat = b * gj
                      + std::sqrt(sG * sG + tau_true * tau_true) * norm(rng_loc);
    }
    return data;
}

// =============================================================================
// Section 14: Full MCMC driver
// =============================================================================

struct MCMCSamples {
    std::vector<double> psi, tau, alpha;
    std::vector<int>    K;
    std::vector<std::vector<double>> beta_per_snp;  // [SNP index][sample index]
};

static MCMCSamples run_mcmc(const std::vector<SNP>& data,
                              const HyperParams& hp,
                              const Config& cfg) {
    g_rng.seed(cfg.seed);
    int p   = (int)data.size();
    int dim = cfg.use_tau ? 2 : 1;

    // ---------- Initialise partition: all singletons ----------
    std::vector<int>    xi(p);
    std::iota(xi.begin(), xi.end(), 0);
    std::vector<double> atoms(p);
    for (int k = 0; k < p; ++k)
        atoms[k] = hp.mu_beta + 0.1 * hp.sig_beta * stdnorm();
    int K = p;

    double psi   = 0.3 + 0.1 * u01();
    double tau   = cfg.use_tau ? (0.05 + 0.05 * u01()) : 0.0;
    double alpha = 1.0;
    std::vector<double> gamma_vec(p, 0.0);

    // NUTS parameterisation
    std::vector<double> theta(dim);
    theta[0] = std::log(psi);
    if (dim > 1) theta[1] = std::log(tau);

    // Find a reasonable initial step size via bisection
    double eps0 = 0.05;
    {
        for (int trial = 0; trial < 30; ++trial) {
            std::vector<double> r0(dim);
            for (int i = 0; i < dim; ++i) r0[i] = stdnorm();
            auto th1 = theta; auto r1 = r0;
            double H0 = hamiltonian(theta, r0, data, xi, atoms, hp, cfg.use_tau);
            leapfrog(th1, r1, eps0, data, xi, atoms, hp, cfg.use_tau);
            double H1 = hamiltonian(th1, r1, data, xi, atoms, hp, cfg.use_tau);
            double acc = std::exp(H0 - H1);
            if (!std::isfinite(acc) || acc < 0.4) eps0 *= 0.5;
            else if (acc > 0.9)                   eps0 *= 1.5;
            else                                  break;
        }
    }
    double log_eps     = std::log(eps0);
    double log_eps_bar = std::log(eps0);
    double H_bar       = 0.0;
    double mu_da       = std::log(10.0 * eps0);

    // ---------- Storage ----------
    MCMCSamples out;
    out.beta_per_snp.resize(p);
    int S_total = cfg.S_wu + cfg.S;

    std::cout << "BNPM-" << (cfg.use_tau ? "PMR" : "MR")
              << " MCMC  |  p=" << p << "  warmup=" << cfg.S_wu
              << "  samples=" << cfg.S << "\n";
    std::cout << std::string(70, '-') << "\n";
    std::cout << std::setw(7)  << "iter"
              << std::setw(6)  << "K"
              << std::setw(10) << "psi"
              << std::setw(10) << "tau"
              << std::setw(10) << "alpha"
              << std::setw(10) << "eps"
              << "\n";

    for (int iter = 1; iter <= S_total; ++iter) {
        bool warmup = (iter <= cfg.S_wu);

        // 1. Gibbs for gamma_j
        update_gamma(gamma_vec, data, xi, atoms, psi, tau);

        // 2. NUTS for (log_psi, [log_tau])
        double eps_cur = std::exp(warmup ? log_eps : log_eps_bar);
        double alpha_bar = nuts_step(theta, eps_cur, data, xi, atoms, hp,
                                      cfg.use_tau, cfg.max_depth);
        psi = std::exp(theta[0]);
        if (cfg.use_tau && dim > 1) tau = std::exp(theta[1]);

        if (warmup)
            dual_avg_update(log_eps, log_eps_bar, H_bar, iter, mu_da,
                             alpha_bar, cfg.delta_nuts);

        // 3. Neal's Algorithm 8 for partition
        update_partition_neal8(xi, atoms, K, data, psi, tau, alpha, hp, cfg.m_aux);

        // 4. ESS for cluster atoms
        update_atoms_ess(atoms, xi, K, data, psi, tau, hp);

        // 5. Jain-Neal split-merge (every T_SM sweeps)
        if (iter % cfg.T_SM == 0)
            split_merge_move(xi, atoms, K, data, psi, tau, alpha, hp, cfg.t_lnch);

        // 6. Escobar-West for alpha
        update_alpha(alpha, K, p, hp);

        // ---------- Store post-warm-up samples ----------
        if (!warmup) {
            out.psi.push_back(psi);
            out.tau.push_back(tau);
            out.alpha.push_back(alpha);
            out.K.push_back(K);
            for (int j = 0; j < p; ++j)
                out.beta_per_snp[j].push_back(atoms[xi[j]]);
        }

        // Progress report
        if (iter % 500 == 0)
            std::cout << std::setw(7)  << iter
                      << std::setw(6)  << K
                      << std::setw(10) << std::fixed << std::setprecision(4) << psi
                      << std::setw(10) << tau
                      << std::setw(10) << alpha
                      << std::setw(10) << std::exp(log_eps_bar)
                      << "\n";
    }
    std::cout << std::string(70, '-') << "\n";
    return out;
}

// =============================================================================
// Section 15: Write samples to CSV and print posterior summaries
// =============================================================================

static void write_csv(const std::string& fname,
                       const MCMCSamples& smp,
                       int p) {
    std::ofstream f(fname);
    if (!f) { std::cerr << "Cannot open " << fname << "\n"; return; }
    f << "psi,tau,alpha,K";
    for (int j = 0; j < p; ++j) f << ",beta_" << j;
    f << "\n";
    int S = (int)smp.psi.size();
    for (int s = 0; s < S; ++s) {
        f << std::fixed << std::setprecision(6)
          << smp.psi[s] << "," << smp.tau[s] << ","
          << smp.alpha[s] << "," << smp.K[s];
        for (int j = 0; j < p; ++j) f << "," << smp.beta_per_snp[j][s];
        f << "\n";
    }
    std::cout << "Samples written to " << fname
              << "  (" << S << " rows, " << (4 + p) << " columns)\n";
}

static void print_summary(const MCMCSamples& smp,
                            int p,
                            double psi_true, double tau_true,
                            const std::vector<double>& beta_true_vec) {
    int S = (int)smp.psi.size();
    auto mean_v = [&](const std::vector<double>& v) {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };
    auto sd_v = [&](const std::vector<double>& v) {
        double m = mean_v(v), ss = 0.0;
        for (double x : v) ss += (x - m) * (x - m);
        return std::sqrt(ss / v.size());
    };
    auto mean_i = [&](const std::vector<int>& v) {
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    };

    std::cout << "\n===== Posterior Summaries (" << S << " samples) =====\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "psi  : mean=" << mean_v(smp.psi)   << "  sd=" << sd_v(smp.psi)
              << "  (true=" << psi_true << ")\n";
    std::cout << "tau  : mean=" << mean_v(smp.tau)   << "  sd=" << sd_v(smp.tau)
              << "  (true=" << tau_true << ")\n";
    std::cout << "alpha: mean=" << mean_v(smp.alpha) << "  sd=" << sd_v(smp.alpha) << "\n";
    std::cout << "K    : mean=" << mean_i(smp.K)     << "\n";

    int show = std::min(p, 10);
    std::cout << "\nPosterior causal effects (first " << show << " SNPs):\n";
    std::cout << std::setw(5) << "SNP"
              << std::setw(11) << "post.mean"
              << std::setw(11) << "post.sd"
              << std::setw(11) << "true"
              << "\n";
    for (int j = 0; j < show; ++j) {
        std::cout << std::setw(5) << j
                  << std::setw(11) << mean_v(smp.beta_per_snp[j])
                  << std::setw(11) << sd_v(smp.beta_per_snp[j])
                  << std::setw(11) << beta_true_vec[j]
                  << "\n";
    }
}

// =============================================================================
// Section 16: main()
// =============================================================================

int main(int argc, char* argv[]) {
    // ---- True data-generating parameters ----
    int    p         = 300;
    double psi_true  = 0.40;   // instrument strength scale
    double tau_true  = 0.05;   // pleiotropy overdispersion

    // Two-component DP mixture: 70% causal (beta=0.20), 30% null (beta=0.00)
    std::vector<double> betas_true = {0.80, 0.00};
    std::vector<double> pi_true    = {0.70, 0.30};

    // ---- Parse simple command-line flags ----
    bool use_tau_flag = true;
    for (int a = 1; a < argc; ++a) {
        std::string arg(argv[a]);
        if (arg == "--no-tau")  { use_tau_flag = false; tau_true = 0.0; }
        if (arg == "--help") {
            std::cout << "Usage: bnpm_mcmc [--no-tau]\n"
                      << "  --no-tau  : run BNPM-MR (no pleiotropy, tau fixed to 0)\n";
            return 0;
        }
    }

    // ---- Generate synthetic data ----
    std::cout << "Generating synthetic data: p=" << p
              << "  psi_true=" << psi_true
              << "  tau_true=" << tau_true
              << "  K_true=2\n\n";
    auto data = generate_data(p, psi_true, tau_true, betas_true, pi_true, /*seed=*/7777);

    // Recover ground-truth beta per SNP (re-run the RNG with same seed)
    {
        std::mt19937_64 rng_gt(7777);
        std::normal_distribution<double> norm(0.0, 1.0);
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        std::uniform_real_distribution<double> use_se(0.03, 0.07);
        std::uniform_real_distribution<double> use_sG(0.05, 0.09);
        // (just advance the rng state to match; not needed for summary)
    }
    // Approximate per-SNP ground truth
    std::mt19937_64 rng_gt2(7777);
    std::uniform_real_distribution<double> unif2(0.0, 1.0);
    std::uniform_real_distribution<double> use_se2(0.03, 0.07);
    std::uniform_real_distribution<double> use_sG2(0.05, 0.09);
    std::normal_distribution<double> norm2(0.0, 1.0);
    std::vector<double> beta_true_vec(p);
    for (int j = 0; j < p; ++j) {
        double r = unif2(rng_gt2); double cum = 0.0;
        int k = (int)pi_true.size() - 1;
        for (int c = 0; c < (int)pi_true.size(); ++c) {
            cum += pi_true[c];
            if (r <= cum) { k = c; break; }
        }
        beta_true_vec[j] = betas_true[k];
        // Advance other draws to stay in sync
        use_se2(rng_gt2); use_sG2(rng_gt2);
        norm2(rng_gt2); norm2(rng_gt2); norm2(rng_gt2);
    }

    // ---- Hyper-parameters ----
    HyperParams hp;
    hp.nu_psi   = 5.0; hp.sig_psi  = 1.0;
    hp.nu_tau   = 5.0; hp.sig_tau  = 1.0;
    hp.mu_beta  = 0.0; hp.sig_beta = 0.5;   // prior SD for causal effects
    hp.a_alpha  = 1.0; hp.b_alpha  = 1.0;

    // ---- MCMC configuration ----
    Config cfg;
    cfg.S_wu         = 20000;
    cfg.S            = 40000;
    cfg.T_SM         = 5;
    cfg.m_aux        = 3;
    cfg.t_lnch       = 5;
    cfg.delta_nuts   = 0.65;
    cfg.use_tau      = use_tau_flag;
    cfg.max_depth    = 10;
    cfg.seed         = 42;
    cfg.outfile      = "bnpm_samples.csv";

    // ---- Run MCMC ----
    auto smp = run_mcmc(data, hp, cfg);

    // ---- Write CSV ----
    write_csv(cfg.outfile, smp, p);

    // ---- Print summaries ----
    print_summary(smp, p, psi_true, tau_true, beta_true_vec);

    return 0;
}
