// mcmc.cpp

#include "bayesmr.h"

void bayesmr_mcmc_noclus(
  double* gamma_chain,
  double* beta_chain,
  double* accept,
  double* loglik,
  double* logprior,
  double* logpost,
  double* data,
  double gamma_p,                             // gamma starting value
  double beta_p,                              // beta starting value
  const double hyper_gammaj_psi2,
  const double rhyper_Gammaj_tau2,
  const double rhyper_gamma_mean,
  const double rhyper_gamma_var,
  const double rhyper_beta_mean,
  const double rhyper_beta_var,
  const double sigma2_beta,                   // beta proposal variance (note sqrt used)
  int totiter,
  int n,
  int p,
  int G,
  int verbose){
  // Sanity
  if (n <= 0 || totiter <= 0) return;

  // Convert input data into vectors once (faster & safer)
  std::vector<double> gammahat_j(n);
  std::vector<double> Gammahat_j(n);
  std::vector<double> sigma2_X(n);
  std::vector<double> sigma2_Y(n);

  // Pre-allocate arrays reused every iteration
  std::vector<double> psi2_j(n);
  std::vector<double> tau2_j(n);
  std::vector<double> a_j(n);
  std::vector<double> v_j(n);

  for (int i = 0; i < n; ++i) {
    gammahat_j[i] = data[i];
    Gammahat_j[i] = data[n + i];
    double sx = data[2 * n + i];
    double sy = data[3 * n + i];
    sigma2_X[i] = sx * sx;
    sigma2_Y[i] = sy * sy;
    psi2_j[i] = sigma2_X[i] + hyper_gammaj_psi2;
    tau2_j[i] = sigma2_Y[i] + rhyper_Gammaj_tau2;
  }

  // RNG
  GetRNGstate();

  // bookkeeping
  double beta_old = beta_p;
  double accept_beta = 0.0;   // track total accept for this chain
  accept[0] = 0.0;

  // local temporaries
  const double hyper_gamma_var_inv = 1.0 / rhyper_gamma_var;
  const double sqrt_rhyper_gamma_var = std::sqrt(rhyper_gamma_var);
  const double sqrt_rhyper_beta_var = std::sqrt(rhyper_beta_var);
  const double sqrt_sigma2_beta = std::sqrt(sigma2_beta);

  // main MCMC loop
  for (int iter = 1; iter <= totiter; ++iter) {
    // -------------------------
    // 1) Gibbs update for gamma
    // -------------------------
    // Compute sums needed for A_beta and B_beta efficiently, avoiding repeated pow()
    const double beta2 = beta_old * beta_old;
    const double beta_inv = (beta_old != 0.0) ? 1.0 / beta_old : 0.0; // used for sum_B_beta_2

    double sum_A_beta_1 = 0.0;
    double sum_A_beta_2 = 0.0;
    double sum_B_beta_1 = 0.0;
    double sum_B_beta_2 = 0.0;

    bool invalid_v = false;
    for (int i = 0; i < n; ++i) {
      a_j[i] = beta2 * hyper_gammaj_psi2 + tau2_j[i];
      const double c = beta_old * hyper_gammaj_psi2;
      const double c2 = c * c;
      v_j[i] = a_j[i] * psi2_j[i] - c2;

      // If determinant non-positive -> mark invalid (posterior undefined); handle gracefully
      if (v_j[i] <= 0.0) {
        invalid_v = true;
        break;
      }

      // accumulate
      const double inv_v = 1.0 / v_j[i];
      sum_A_beta_1 += tau2_j[i] * inv_v;
      sum_A_beta_2 += sigma2_X[i] * inv_v;
      sum_B_beta_1 += gammahat_j[i] * tau2_j[i] * inv_v;
      // guard against division by zero on beta_old
      if (beta_old != 0.0) {
        sum_B_beta_2 += (Gammahat_j[i] * sigma2_X[i]) * (inv_v * beta_inv);
      } else {
        // if beta_old == 0, the expression becomes infinite/undefined; set flag
        invalid_v = true;
        break;
      }
    }

    // If matrix was singular or beta==0 caused issues, reject/handle:
    // For simplicity: if invalid, we skip the Gibbs update and keep previous gamma
    double gamma_new = gamma_p; // fallback (starting gamma)
    if (!invalid_v) {
      const double A_beta = sum_A_beta_1 + beta2 * sum_A_beta_2 + hyper_gamma_var_inv;
      const double B_beta = sum_B_beta_1 + beta2 * sum_B_beta_2 + (rhyper_gamma_mean * hyper_gamma_var_inv);
      // Draw gamma ~ N(mean = B/A, var = 1/A)
      const double mean_gamma = B_beta / A_beta;
      const double sd_gamma = std::sqrt(1.0 / A_beta);
      gamma_new = R::rnorm(mean_gamma, sd_gamma);
    } else {
      // keep previous gamma (or you might want to set to prior mean)
      // We choose to keep previous gamma (gamma_p is starting value; or you can track last sampled)
      // For clarity, set gamma_new to last stored gamma in chain if available:
      if ( (iter - 2) >= 0 ) {
          gamma_new = gamma_chain[iter - 2];
      } else {
          gamma_new = gamma_p;
      }
    }

    gamma_chain[iter - 1] = gamma_new;

    // -------------------------------------------
    // 2) Metropolis-Hastings update for beta
    // -------------------------------------------
    // Propose new beta (random-walk Normal)
    double beta_prop = R::rnorm(beta_old, sqrt_sigma2_beta);

    // compute log-posterior for proposed and current betas
    double lpost_prop = 0.0;
    double lpost_old = 0.0;

    logpost_beta(&lpost_prop, beta_prop, gamma_new,
      rhyper_beta_mean, rhyper_beta_var, hyper_gammaj_psi2, rhyper_Gammaj_tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    logpost_beta(&lpost_old, beta_old, gamma_new,
      rhyper_beta_mean, rhyper_beta_var, hyper_gammaj_psi2, rhyper_Gammaj_tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    // compute acceptance probability (handle -inf or NaN)
    double prob = 0.0;
    if (!std::isfinite(lpost_prop)) {
      prob = 0.0;
    } else if (!std::isfinite(lpost_old)) {
      // if old is -inf but new is finite, accept
      prob = 1.0;
    } else {
      const double diff = lpost_prop - lpost_old;
      // guard against overflow
      if (diff >= 0.0) prob = 1.0;
      else prob = std::exp(diff);
    }

    double u = R::runif(0.0, 1.0);
    int accepted_beta = (u < prob) ? 1 : 0;
    if (accepted_beta) {
      beta_old = beta_prop;
      accept_beta += 1.0;
    }
    accept[0] += accepted_beta;   // cumulative acceptances (user originally used accept[0])

    // store beta
    beta_chain[iter - 1] = beta_old;

    // ---------------------------------------------------
    // 3) compute log-densities: logprior, loglik, logpost
    // ---------------------------------------------------
    // logprior: normal for gamma and beta (log-scale)
    const double logprior_gamma = R::dnorm(gamma_new, rhyper_gamma_mean, sqrt_rhyper_gamma_var, 1);
    const double logprior_beta  = R::dnorm(beta_old, rhyper_beta_mean, sqrt_rhyper_beta_var, 1);
    logprior[iter - 1] = logprior_gamma + logprior_beta;

    // loglik using optimized bayesmr_logLik
    extern double bayesmr_logLik(const double beta, const double gamma,
      const double psi2, const double tau2, int n,
      const double* gammahat_j, const double* Gammahat_j,
      const double* sigma2_X, const double* sigma2_Y);

    loglik[iter - 1] = bayesmr_logLik(beta_old, gamma_new,
      hyper_gammaj_psi2, rhyper_Gammaj_tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    // logpost
    logpost[iter - 1] = loglik[iter - 1] + logprior[iter - 1];

    // optional progress printing
    if ((iter % 500) == 0 && verbose) {
      double accept_rate = accept_beta / static_cast<double>(G*iter);
      REprintf("   iteration %d/%d ==> acceptance beta: %1.4f\n", iter, totiter, accept_rate);
    }

    // interruption point (allow user to interrupt R)
    R_CheckUserInterrupt();
  } // end iterations

  PutRNGstate();

  // done - vectors auto-clean
}
