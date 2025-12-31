// mcmc_het.cpp

#include "bayesmr.h"

void bayesmr_mcmc_noclus_het(
  double* gamma_chain,
  double* beta_chain,
  double* psi_chain,
  double* tau_chain,
  double* accept,
  double* loglik,
  double* logprior,
  double* logpost,
  double* data,
  double gamma_p,                             // gamma starting value
  double beta_p,                              // beta starting value
  double psi_p,                               // psi starting value
  double tau_p,                               // tau starting value
  const double rhyper_alpha_psi,
  const double rhyper_nu_psi,
  const double rhyper_alpha_tau,
  const double rhyper_nu_tau,
  const double rhyper_gamma_mean,
  const double rhyper_gamma_var,
  const double rhyper_beta_mean,
  const double rhyper_beta_var,
  const double sigma2_beta,                   // beta proposal variance (note sqrt used)
  const double C_psi,                         // psi proposal parameter
  const double C_tau,                         // tau proposal parameter
  int totiter,
  int n,
  int p,
  int G,
  int verbose){
  if (n <= 0 || totiter <= 0) return;

  // Convert input data into vectors once (faster & safer)
  std::vector<double> gammahat_j(n);
  std::vector<double> Gammahat_j(n);
  std::vector<double> sigma2_X(n);
  std::vector<double> sigma2_Y(n);

  for (int i = 0; i < n; ++i) {
    gammahat_j[i] = data[i];
    Gammahat_j[i] = data[n + i];
    double sx = data[2 * n + i];
    double sy = data[3 * n + i];
    sigma2_X[i] = sx * sx;
    sigma2_Y[i] = sy * sy;
  }

  // RNG
  GetRNGstate();

  double gamma_old = gamma_p;  // [unnecessary]
  double beta_old = beta_p;
  double psi_old = psi_p;
  double tau_old = tau_p;

  double accept_beta = 0.0;    // track total accept for beta chain
  double accept_psi = 0.0;     // track total accept for psi chain
  double accept_tau = 0.0;     // track total accept for tau chain

  double beta2 = beta_old * beta_old;
  double psi2 = psi_old * psi_old;
  double tau2 = tau_old * tau_old;

  // local temporaries
  const double rhyper_gamma_var_inv = 1.0 / rhyper_gamma_var;
  const double sqrt_rhyper_gamma_var = std::sqrt(rhyper_gamma_var);
  const double sqrt_rhyper_beta_var = std::sqrt(rhyper_beta_var);
  const double sqrt_sigma2_beta = std::sqrt(sigma2_beta);

  // Fixed proposal standard deviations (log scale)
  const double SD_LOG_ETA   = 0.15;   // ~15% multiplicative moves
  const double SD_LOG_OMEGA = 0.15;

  // main MCMC loop
  for (int iter = 1; iter <= totiter; ++iter) {
    // -------------------------
    // 1) Gibbs update for gamma
    // -------------------------
    double sum_A_beta_1 = 0.0, sum_A_beta_2 = 0.0, sum_B_beta_1 = 0.0, sum_B_beta_2 = 0.0;

    for (int i = 0; i < n; ++i) {
      double psi2_j = sigma2_X[i] + psi2;
      double tau2_j = sigma2_Y[i] + tau2;

      double a_j = beta2 * psi2 + tau2_j;
      double c = beta_old * psi2;
      double v_j = a_j * psi2_j - c * c;
      if (v_j <= 0){ return; }  // exit if variance is not positive

      // accumulate
      const double inv_v_j = 1.0 / v_j;
      sum_A_beta_1 += tau2_j * inv_v_j;
      sum_A_beta_2 += sigma2_X[i] * inv_v_j;
      sum_B_beta_1 += gammahat_j[i] * tau2_j * inv_v_j;
      sum_B_beta_2 += (Gammahat_j[i] * sigma2_X[i]) * inv_v_j;
    }

    const double A_beta = sum_A_beta_1 + beta2 * sum_A_beta_2 + rhyper_gamma_var_inv;
    const double B_beta = sum_B_beta_1 + beta_old * sum_B_beta_2 + (rhyper_gamma_mean * rhyper_gamma_var_inv);
    // Draw gamma ~ N(mean = B/A, var = 1/A)
    const double mean_gamma = B_beta / A_beta;
    const double sd_gamma = std::sqrt(1.0 / A_beta);
    gamma_old = R::rnorm(mean_gamma, sd_gamma);

    gamma_chain[iter - 1] = gamma_old;

    // --------------------------------------
    // 2) Metropolis-Hastings update for beta
    // --------------------------------------
    // Propose new beta (random-walk Normal)
    double beta_prop = R::rnorm(beta_old, sqrt_sigma2_beta);

    // compute log-posterior for proposed and current betas
    double lpost_prop = 0.0, lpost_curr = 0.0;

    logpost_beta(&lpost_prop, beta_prop, gamma_old,
      rhyper_beta_mean, rhyper_beta_var, psi2, tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    logpost_beta(&lpost_curr, beta_old, gamma_old,
      rhyper_beta_mean, rhyper_beta_var, psi2, tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    // compute acceptance probability (handle -inf or NaN)
    double prob = 0.0;
    if (!std::isfinite(lpost_prop)) {
      prob = 0.0;
    } else if (!std::isfinite(lpost_curr)) {
      // if old is -inf but new is finite, accept
      prob = 1.0;
    } else {
      const double diff = lpost_prop - lpost_curr;
      // guard against overflow
      if (diff >= 0.0) prob = 1.0;
      else prob = std::exp(diff);
    }

    double u = R::runif(0.0, 1.0);
    if (u < prob) {
      beta_old = beta_prop;
      beta2 = beta_old * beta_old;
      accept_beta++;
    }

    // store beta
    beta_chain[iter - 1] = beta_old;

    // -------------------------------------
    // 3+4) Fixed MH update in (eta, omega)
    //      (no adaptation)
    // -------------------------------------

    const double beta_abs = std::fabs(beta_old);

    // current transformed parameters
    double eta_old   = beta_abs * psi_old;
    double omega_old = std::sqrt(tau2 + eta_old * eta_old);

    // log-space
    double zeta_old = std::log(eta_old + 1e-12);
    double zow_old  = std::log(omega_old);

    // fixed independent RW proposals
    double zeta_prop = R::rnorm(zeta_old, SD_LOG_ETA);
    double zow_prop  = R::rnorm(zow_old,  SD_LOG_OMEGA);

    // back-transform
    double eta_prop   = std::exp(zeta_prop);
    double omega_prop = std::exp(zow_prop);

    // placeholders
    bool accept_etaomega = false;
    double psi_prop = psi_old;
    double tau_prop = tau_old;

    // feasibility check
    double tmp = omega_prop * omega_prop - eta_prop * eta_prop;

    if (tmp > 0.0 && beta_abs > 0.0) {

      psi_prop = eta_prop / beta_abs;
      tau_prop = std::sqrt(tmp);

      // ---- proposed log-posterior ----
      double lp_prop =
        bayesmr_logLik(beta_old, gamma_old,
                       psi_prop * psi_prop,
                       tau_prop * tau_prop,
                       n,
                       gammahat_j.data(),
                       Gammahat_j.data(),
                       sigma2_X.data(),
                       sigma2_Y.data());

      lp_prop += dhalft_scalar(psi_prop,
                               rhyper_alpha_psi,
                               rhyper_nu_psi,
                               true);

      lp_prop += dhalft_scalar(tau_prop,
                               rhyper_alpha_tau,
                               rhyper_nu_tau,
                               true);

      // ---- current log-posterior ----
      double lp_curr =
        bayesmr_logLik(beta_old, gamma_old,
                       psi2, tau2,
                       n,
                       gammahat_j.data(),
                       Gammahat_j.data(),
                       sigma2_X.data(),
                       sigma2_Y.data());

      lp_curr += dhalft_scalar(psi_old,
                               rhyper_alpha_psi,
                               rhyper_nu_psi,
                               true);

      lp_curr += dhalft_scalar(tau_old,
                               rhyper_alpha_tau,
                               rhyper_nu_tau,
                               true);

      // ---- Jacobian correction ----
      // |J| = omega / (|beta| * tau)
      // beta cancels in MH ratio
      double logJ_prop = std::log(omega_prop) - std::log(tau_prop);
      double logJ_curr = std::log(omega_old)  - std::log(tau_old);

      double diff = (lp_prop + logJ_prop) - (lp_curr + logJ_curr);
      double prob = (diff >= 0.0 ? 1.0 : std::exp(diff));

      if (R::runif(0.0, 1.0) < prob) {
        accept_etaomega = true;
      }
    }

    // ---- commit if accepted ----
    if (accept_etaomega) {
      psi_old = psi_prop;
      tau_old = tau_prop;
      psi2 = psi_old * psi_old;
      tau2 = tau_old * tau_old;
      accept_psi++;
      accept_tau++;
    }

    // ---- store ----
    psi_chain[iter - 1] = psi_old;
    tau_chain[iter - 1] = tau_old;

    // ---------------------------------------------------
    // 5) compute log-densities: logprior, loglik, logpost
    // ---------------------------------------------------
    // logprior: normal for gamma and beta (log-scale)
    const double logprior_gamma = R::dnorm(gamma_old, rhyper_gamma_mean, sqrt_rhyper_gamma_var, 1);
    const double logprior_beta  = R::dnorm(beta_old, rhyper_beta_mean, sqrt_rhyper_beta_var, 1);
    const double logprior_psi   = dhalft(std::vector<double>{psi_old},
      std::vector<double>{rhyper_alpha_psi}, std::vector<double>{rhyper_nu_psi}, true)[0];
    const double logprior_tau   = dhalft(std::vector<double>{tau_old},
      std::vector<double>{rhyper_alpha_tau}, std::vector<double>{rhyper_nu_tau}, true)[0];
    logprior[iter - 1] = logprior_gamma + logprior_beta + logprior_psi + logprior_tau;

    loglik[iter - 1] = bayesmr_logLik(beta_old, gamma_old, psi2, tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    // logpost
    logpost[iter - 1] = loglik[iter - 1] + logprior[iter - 1];

    // optional progress printing
    if (verbose && (iter % 500 == 0)){
      REprintf("iter %d/%d ==> acc(beta) = %.3f - acc(psi) = %.3f - acc(tau) = %.3f\n",
        iter, totiter,
        accept_beta/iter,
        accept_psi/iter,
        accept_tau/iter);
    }

    // interruption point (allow user to interrupt R)
    R_CheckUserInterrupt();
  } // end iterations

  accept[0] = accept_beta / totiter;
  accept[1] = accept_psi  / totiter;
  accept[2] = accept_tau  / totiter;

  PutRNGstate();

  // done - vectors auto-clean
}
