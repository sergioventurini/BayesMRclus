// mcmc_mix_het_2.cpp (includes psi/tau reparameterization)

// #include "bayesmr.h"

// void bayesmr_mcmc_mix_het(
//   double* gamma_chain,
//   double* beta_chain,
//   int* xi_chain,
//   double* alpha_chain,
//   double* psi_chain,
//   double* tau_chain,
//   double* accept,
//   double* loglik,
//   double* logprior,
//   double* logpost,
//   double* data,
//   double gamma_p,                             // gamma starting value
//   double* beta_p,                             // beta starting values
//   int* xi_p,                                  // xi starting values
//   double alpha_p,                             // alpha starting value
//   double psi_p,                               // psi starting value
//   double tau_p,                               // tau starting value
//   double K_p,                                 // K (# clusters) starting value
//   const double rhyper_gamma_mean,
//   const double rhyper_gamma_var,
//   const double rhyper_beta_mean,
//   const double rhyper_beta_var,
//   const double rhyper_psi_alpha,
//   const double rhyper_psi_nu,
//   const double rhyper_tau_alpha,
//   const double rhyper_tau_nu,
//   const double rhyper_alpha_a,
//   const double rhyper_alpha_b,
//   const double C_beta,                        // beta proposal sd
//   const double C_logpsi,                      // log-psi proposal sd
//   const double C_logtau,                      // log-tau proposal sd
//   int m_beta,                                 // # of proposed beta values
//   int totiter,
//   int n,
//   int verbose){
//   if (n <= 0 || totiter <= 0) return;

//   // Convert input data into vectors once (faster & safer)
//   std::vector<double> gammahat_j(n);
//   std::vector<double> Gammahat_j(n);
//   std::vector<double> sigma2_X(n);
//   std::vector<double> sigma2_Y(n);

//   for (size_t i = 0; i < n; ++i) {
//     gammahat_j[i] = data[i];
//     Gammahat_j[i] = data[n + i];
//     double sx = data[2 * n + i];
//     double sy = data[3 * n + i];
//     sigma2_X[i] = sx * sx;
//     sigma2_Y[i] = sy * sy;
//   }

//   // RNG
//   GetRNGstate();

//   size_t K_old = K_p;
//   double gamma_old = gamma_p;     // [unnecessary]
//   std::vector<double> beta_star_old(K_old);
//   for (size_t j = 0; j < K_old; ++j) {
//     beta_star_old[j] = beta_p[j];
//   }
//   std::vector<int> xi_old(n);
//   for (size_t j = 0; j < n; ++j) {
//     xi_old[j] = xi_p[j] - 1;      // uses 0-based indexing
//   }
//   double alpha_old = alpha_p;     // [unnecessary]
//   double accept_beta = 0.0;       // track total accept for beta chain
//   double psi_old = psi_p;
//   double accept_psi = 0.0;        // track total accept for psi chain
//   double tau_old = tau_p;
//   double accept_tau = 0.0;        // track total accept for tau chain
//   double total_beta_proposals = 0.0;
//   double total_psitau_proposals = 0.0;

//   // local temporaries
//   const double sqrt_rhyper_gamma_var = std::sqrt(rhyper_gamma_var);
//   const double sqrt_rhyper_beta_var = std::sqrt(rhyper_beta_var);

//   // main MCMC loop
//   for (int iter = 1; iter <= totiter; ++iter) {
//     // ------------------
//     // 1) Update of gamma
//     // ------------------
//     gamma_old = full_cond_gamma(rhyper_gamma_mean, rhyper_gamma_var,
//       gammahat_j, Gammahat_j, sigma2_X, sigma2_Y, 
//       psi_old, tau_old, beta_star_old, xi_old);

//     gamma_chain[iter - 1] = gamma_old;

//     // ---------------
//     // 2) Update of xi
//     // ---------------
//     K_old = beta_star_old.size();  // allocated clusters
//     std::vector<std::vector<double>> weights(n, std::vector<double>(K_old + m_beta));

//     std::vector<int> cluster_counts(K_old, 0);
//     for (size_t h = 0; h < n; ++h) cluster_counts[xi_old[h]]++;

//     for (size_t k = 0; k < K_old; ++k) {
//       for (size_t j = 0; j < n; ++j) {
//         double ghat = gammahat_j[j];
//         double Ghat = Gammahat_j[j];
//         double s2x  = sigma2_X[j];
//         double s2y  = sigma2_Y[j];
//         double loglik_jk = bayesmr_logLik(
//           beta_star_old[k], gamma_old, psi_old*psi_old,
//           tau_old*tau_old, 1, &ghat, &Ghat, &s2x, &s2y
//         );

//         // int n_jk = count_excluding(xi_old, k, j);
//         int n_jk = cluster_counts[k] - (xi_old[j] == k ? 1 : 0);
//         weights[j][k] = std::log(static_cast<double>(n_jk)) + loglik_jk;
//       }
//     }

//     // Draw m new betas from the baseline measure
//     std::vector<double> beta_new(m_beta);
//     for (int m = 0; m < m_beta; ++m) {
//       beta_new[m] = R::rnorm(rhyper_beta_mean, sqrt_rhyper_beta_var);
//     }

//     for (int m = 0; m < m_beta; ++m) {
//       for (size_t j = 0; j < n; ++j) {
//         double ghat = gammahat_j[j];
//         double Ghat = Gammahat_j[j];
//         double s2x  = sigma2_X[j];
//         double s2y  = sigma2_Y[j];
//         double loglik_jk = bayesmr_logLik(
//           beta_new[m], gamma_old, psi_old*psi_old,
//           tau_old*tau_old, 1, &ghat, &Ghat, &s2x, &s2y
//         );

//         weights[j][K_old + m] = std::log(alpha_old / m_beta) + loglik_jk;
//       }
//     }

//     // Normalize and update xi by sampling from 1:(K + m) with probabilities given by the j-th row of probs
//     for (size_t j = 0; j < n; ++j) {
//       std::vector<double> probs = normalize_weights(weights[j]);
//       xi_old[j] = sample_discrete(probs);
//     }

//     // Update beta_star
//     std::vector<double> betas = beta_star_old;
//     betas.insert(betas.end(), beta_new.begin(), beta_new.end());

//     std::set<int> unique_xi(xi_old.begin(), xi_old.end());
//     std::vector<double> beta_star_new;
//     for (int k : unique_xi) {
//       beta_star_new.push_back(betas[k]);
//     }
//     beta_star_old = beta_star_new;

//     relabel_xi(xi_old);
//     K_old = beta_star_old.size();

//     for (size_t j = 0; j < n; ++j) {
//       xi_chain[(iter - 1)*n + j] = xi_old[j] + 1;   // return 1-based indexes
//     }

//     // --------------------------------------------------
//     // 3) Metropolis-Hastings update of each beta_star[k]
//     // --------------------------------------------------
//     for (size_t k = 0; k < K_old; ++k) {
//       double beta_k = beta_star_old[k];
//       double beta_k_prop = R::rnorm(beta_k, C_beta);
//       total_beta_proposals++;

//       auto gamma_hat_k = subset(gammahat_j, xi_old, k);
//       auto Gamma_hat_k = subset(Gammahat_j, xi_old, k);
//       auto sigma2_X_k = subset(sigma2_X, xi_old, k);
//       auto sigma2_Y_k = subset(sigma2_Y, xi_old, k);

//       double r_beta_k = full_cond_beta(beta_k_prop, rhyper_beta_mean, rhyper_beta_var,
//                           gamma_hat_k, Gamma_hat_k, sigma2_X_k, sigma2_Y_k, gamma_old,
//                           psi_old, tau_old) -
//                         full_cond_beta(beta_k, rhyper_beta_mean, rhyper_beta_var,
//                           gamma_hat_k, Gamma_hat_k, sigma2_X_k, sigma2_Y_k, gamma_old,
//                           psi_old, tau_old);

//       double log_u = std::log(unif_rand());
//       if (log_u < r_beta_k) {
//         beta_star_old[k] = beta_k_prop;
//         accept_beta++;
//       }
//     }

//     for (size_t j = 0; j < n; ++j) {
//       if (j < K_old) {
//         beta_chain[(iter - 1)*n + j] = beta_star_old[j];
//       }
//       else {
//         beta_chain[(iter - 1)*n + j] = NA_REAL;
//       }
//     }

//     // --------------------------------------------------
//     // Compute curvature proxy: beta_tilde^2
//     // --------------------------------------------------
//     std::vector<int> cluster_sizes(K_old, 0);
//     for (size_t j = 0; j < n; ++j) {
//       cluster_sizes[xi_old[j]]++;
//     }

//     double beta_tilde_sq = 0.0;
//     for (size_t k = 0; k < K_old; ++k) {
//       beta_tilde_sq += cluster_sizes[k] * beta_star_old[k] * beta_star_old[k];
//     }
//     beta_tilde_sq /= static_cast<double>(n);

//     // safeguard against numerical degeneracy
//     if (beta_tilde_sq < 1e-12) beta_tilde_sq = 1e-12;

//     // ---------------------------------------------------
//     // 4-5) Joint MH update of (psi, tau) via (eta, omega)
//     // ---------------------------------------------------
//     auto beta_long_old = remap_vec(beta_star_old, xi_old);

//     // current eta and omega
//     double eta_old = psi_old * std::sqrt(beta_tilde_sq);
//     double omega_old = std::sqrt(tau_old * tau_old + eta_old * eta_old);

//     // log scale
//     double z_eta_old   = std::log(eta_old);
//     double z_omega_old = std::log(omega_old);

//     // propose random walk
//     double z_eta_prop   = R::rnorm(z_eta_old, C_logpsi);   // reuse tuning
//     double z_omega_prop = R::rnorm(z_omega_old, C_logtau); // reuse tuning

//     double eta_prop   = std::exp(z_eta_prop);
//     double omega_prop = std::exp(z_omega_prop);

//     // constraint: omega >= eta
//     if (omega_prop > eta_prop) {
//       total_psitau_proposals++;

//       // back-transform
//       double psi_prop = eta_prop / std::sqrt(beta_tilde_sq);
//       double tau_sq_prop = omega_prop * omega_prop - eta_prop * eta_prop;
//       double tau_prop = std::sqrt(tau_sq_prop);

//       // log posterior at proposed
//       double logpost_prop =
//         full_cond_psitau(beta_long_old,
//                          rhyper_psi_alpha, rhyper_psi_nu,
//                          rhyper_tau_alpha, rhyper_tau_nu,
//                          gammahat_j, Gammahat_j, sigma2_X, sigma2_Y,
//                          gamma_old, psi_prop, tau_prop);

//       // log posterior at current
//       double logpost_old =
//         full_cond_psitau(beta_long_old,
//                          rhyper_psi_alpha, rhyper_psi_nu,
//                          rhyper_tau_alpha, rhyper_tau_nu,
//                          gammahat_j, Gammahat_j, sigma2_X, sigma2_Y,
//                          gamma_old, psi_old, tau_old);

//       // Jacobian term:
//       // log J = log(psi) + 2*log(omega) - log(tau) = z_eta + 2*z_omega -log(tau)
//       double logJ_prop = z_eta_prop + 2.0*z_omega_prop - std::log(tau_prop);
//       double logJ_old  = z_eta_old  + 2.0*z_omega_old  - std::log(tau_old);

//       double r = (logpost_prop - logpost_old) + (logJ_prop - logJ_old);

//       if (std::log(unif_rand()) < r) {
//         psi_old = psi_prop;
//         tau_old = tau_prop;
//         accept_psi++;   // track jointly
//         accept_tau++;
//       }
//     }

//     psi_chain[iter - 1] = psi_old;
//     tau_chain[iter - 1] = tau_old;

//     // ------------------
//     // 6) Update of alpha
//     // ------------------
//     double eta_dp = R::rbeta(alpha_old + 1, n);
//     while (rhyper_alpha_b <= std::log(eta_dp)) {
//       eta_dp = R::rbeta(alpha_old + 1, n);
//     }

//     double gamma1 = R::rgamma(rhyper_alpha_a + K_old, 1.0 / (rhyper_alpha_b - std::log(eta_dp)));
//     double gamma2 = R::rgamma(rhyper_alpha_a + K_old - 1, 1.0 / (rhyper_alpha_b - std::log(eta_dp)));

//     double prob1 = rhyper_alpha_a + K_old - 1;
//     double prob2 = n*(rhyper_alpha_b - std::log(eta_dp));
//     double total_prob = prob1 + prob2;

//     alpha_old = (unif_rand() < prob1 / total_prob) ? gamma1 : gamma2;

//     alpha_chain[iter - 1] = alpha_old;

//     // ---------------------------------------------------
//     // 7) compute log-densities: logprior, loglik, logpost
//     // ---------------------------------------------------
//     // logprior
//   //   const double logprior_gamma = R::dnorm(gamma_old, rhyper_gamma_mean, sqrt_rhyper_gamma_var, 1);
//   //   const double logprior_beta  = R::dnorm(beta_star_old, rhyper_beta_mean, sqrt_rhyper_beta_var, 1);
//   //   logprior[iter - 1] = logprior_gamma + logprior_beta;
//     logprior[iter - 1] = 0;

//   //   loglik[iter - 1] = bayesmr_logLik(beta_star_old, gamma_old, psi_old*psi_old, tau_old*tau_old, n,
//   //     gammahat_j.data(), Gammahat_j.data(), sigma2_X.data(), sigma2_Y.data());
//     loglik[iter - 1] = 0;

//     // logpost
//     logpost[iter - 1] = loglik[iter - 1] + logprior[iter - 1];

//     // optional progress printing
//     if ((iter % 500) == 0 && verbose) {
//       REprintf("   iteration %d/%d ==> acc. beta: %1.4f - acc. psi: %1.4f - acc. tau: %1.4f\n",
//         iter, totiter, accept_beta / total_beta_proposals, accept_psi / total_psitau_proposals,
//         accept_tau / total_psitau_proposals);
//     }

//     // interruption point (allow user to interrupt R)
//     R_CheckUserInterrupt();
//   } // end iterations

//   accept[0] = accept_beta / total_beta_proposals;
//   accept[1] = accept_psi  / total_psitau_proposals;
//   accept[2] = accept_tau  / total_psitau_proposals;

//   PutRNGstate();

//   // done - vectors auto-clean
// }
