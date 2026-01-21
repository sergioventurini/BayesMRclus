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
  const double rhyper_gammaj_psi2,
  const double rhyper_Gammaj_tau2,
  const double rhyper_gamma_mean,
  const double rhyper_gamma_var,
  const double rhyper_beta_mean,
  const double rhyper_beta_var,
  const double C_beta,                        // beta proposal sd
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
    psi2_j[i] = sigma2_X[i] + rhyper_gammaj_psi2;
    tau2_j[i] = sigma2_Y[i] + rhyper_Gammaj_tau2;
  }

  // RNG
  GetRNGstate();

  double gamma_old = gamma_p;  // [unnecessary]
  double beta_old = beta_p;
  double beta2_old = beta_old * beta_old;
  double accept_beta = 0.0;    // track total accept for beta chain

  // local temporaries
  const double hyper_gamma_var_inv = 1.0 / rhyper_gamma_var;
  const double sqrt_rhyper_gamma_var = std::sqrt(rhyper_gamma_var);
  const double sqrt_rhyper_beta_var = std::sqrt(rhyper_beta_var);

  // main MCMC loop
  for (int iter = 1; iter <= totiter; ++iter) {
    // -------------------------
    // 1) Gibbs update for gamma
    // -------------------------
    double sum_A_beta_1 = 0.0;
    double sum_A_beta_2 = 0.0;
    double sum_B_beta_1 = 0.0;
    double sum_B_beta_2 = 0.0;

    for (int i = 0; i < n; ++i) {
      a_j[i] = beta2_old * rhyper_gammaj_psi2 + tau2_j[i];
      const double c = beta_old * rhyper_gammaj_psi2;
      const double c2 = c * c;
      v_j[i] = a_j[i] * psi2_j[i] - c2;

      // accumulate
      const double inv_v = 1.0 / v_j[i];
      sum_A_beta_1 += tau2_j[i] * inv_v;
      sum_A_beta_2 += sigma2_X[i] * inv_v;
      sum_B_beta_1 += gammahat_j[i] * tau2_j[i] * inv_v;
      sum_B_beta_2 += (Gammahat_j[i] * sigma2_X[i]) * inv_v;
    }

    const double A_beta = sum_A_beta_1 + beta2_old * sum_A_beta_2 + hyper_gamma_var_inv;
    const double B_beta = sum_B_beta_1 + beta_old * sum_B_beta_2 + (rhyper_gamma_mean * hyper_gamma_var_inv);
    // Draw gamma ~ N(mean = B/A, var = 1/A)
    const double mean_gamma = B_beta / A_beta;
    const double sd_gamma = std::sqrt(1.0 / A_beta);
    gamma_old = R::rnorm(mean_gamma, sd_gamma);

    gamma_chain[iter - 1] = gamma_old;

    // --------------------------------------
    // 2) Metropolis-Hastings update for beta
    // --------------------------------------
    // Propose new beta (random-walk Normal)
    double beta_prop = R::rnorm(beta_old, C_beta);

    // compute log-posterior for proposed and current betas
    double lpost_prop = 0.0;
    double lpost_old = 0.0;

    logpost_beta(&lpost_prop, beta_prop, gamma_old,
      rhyper_beta_mean, rhyper_beta_var, rhyper_gammaj_psi2, rhyper_Gammaj_tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    logpost_beta(&lpost_old, beta_old, gamma_old,
      rhyper_beta_mean, rhyper_beta_var, rhyper_gammaj_psi2, rhyper_Gammaj_tau2, n,
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
    if (u < prob) {
      beta_old = beta_prop;
      beta2_old = beta_old * beta_old;
      accept_beta++;
    }

    // store beta
    beta_chain[iter - 1] = beta_old;

    // ---------------------------------------------------
    // 3) compute log-densities: logprior, loglik, logpost
    // ---------------------------------------------------
    // logprior: normal for gamma and beta (log-scale)
    const double logprior_gamma = R::dnorm(gamma_old, rhyper_gamma_mean, sqrt_rhyper_gamma_var, 1);
    const double logprior_beta  = R::dnorm(beta_old, rhyper_beta_mean, sqrt_rhyper_beta_var, 1);
    logprior[iter - 1] = logprior_gamma + logprior_beta;

    loglik[iter - 1] = bayesmr_logLik(beta_old, gamma_old, rhyper_gammaj_psi2, rhyper_Gammaj_tau2, n,
      gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());

    // logpost
    logpost[iter - 1] = loglik[iter - 1] + logprior[iter - 1];

    // optional progress printing
    if ((iter % 500) == 0 && verbose) {
      REprintf("   iteration %d/%d ==> acceptance beta: %1.4f\n", iter, totiter, accept_beta / iter);
    }

    // interruption point (allow user to interrupt R)
    R_CheckUserInterrupt();
  } // end iterations

  accept[0] = accept_beta / totiter;

  PutRNGstate();

  // done - vectors auto-clean
}
