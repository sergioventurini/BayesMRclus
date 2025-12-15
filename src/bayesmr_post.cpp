// bayesmr_post.cpp

#include "bayesmr.h"

// full conditional distribution of the beta parameter
void logpost_beta(double* lpost,
  const double beta, const double gamma,
  const double mu_beta, const double sigma2_beta,
  const double psi2, const double tau2, int n,
  const double* gammahat_j, const double* Gammahat_j,
  const double* sigma2_X, const double* sigma2_Y){
  // Convert raw pointers into vectors once
  std::vector<double> X(gammahat_j, gammahat_j + n);
  std::vector<double> Y(Gammahat_j,  Gammahat_j + n);

  // Precompute beta-related constants
  const double beta2 = beta * beta;
  const double beta_gamma = beta * gamma;
  const double beta_psi2 = beta * psi2;

  // Allocate mean vectors
  std::vector<double> mu_x(n, gamma);
  std::vector<double> mu_y(n, beta_gamma);

  // Allocate covariance vectors
  std::vector<double> sigma_xx(n);
  std::vector<double> sigma_yy(n);
  std::vector<double> sigma_xy(n, beta_psi2);

  // Fill variances per observation
  for (int j = 0; j < n; ++j) {
    sigma_xx[j] = sigma2_X[j] + psi2;
    sigma_yy[j] = beta2 * psi2 + sigma2_Y[j] + tau2;
  }

  // Compute log-density for all observations in ONE call
  std::vector<double> logdens = dbivnorm_cpp(
    X, Y,           // data
    mu_x, mu_y,     // means
    sigma_xx,       // Var(gammahat_j)
    sigma_yy,       // Var(Gammahat_j)
    sigma_xy,       // Cov
    true            // log-scale
  );

  // Sum likelihood
  double logLik = std::accumulate(logdens.begin(), logdens.end(), 0.0);

  // Gaussian prior for beta: -(beta - mu)^2/(2*sigma2)
  const double diff = beta - mu_beta;
  const double logprior = -0.5 * (diff * diff) / sigma2_beta;

  // Final result: log posterior
  *lpost = logLik + logprior;
}

// log-likelihood
// [code optimized to avoid loop-based calculations]
double bayesmr_logLik(const double beta, const double gamma,
  const double psi2, const double tau2, int n,
  const double* gammahat_j, const double* Gammahat_j,
  const double* sigma2_X, const double* sigma2_Y){
  // Convert inputs to vectors once
  std::vector<double> X(gammahat_j, gammahat_j + n);
  std::vector<double> Y(Gammahat_j, Gammahat_j + n);

  // Precompute everything we can
  const double beta2 = beta * beta;
  const double beta_gamma = beta * gamma;
  const double beta_psi2 = beta * psi2;

  // Allocate per-observation parameters
  std::vector<double> mu_x(n, gamma);
  std::vector<double> mu_y(n, beta_gamma);
  std::vector<double> sigma_xx(n);
  std::vector<double> sigma_yy(n);
  std::vector<double> sigma_xy(n, beta_psi2);

  // Fill variance matrix elements
  for (int j = 0; j < n; ++j) {
    sigma_xx[j] = sigma2_X[j] + psi2;
    sigma_yy[j] = beta2 * psi2 + sigma2_Y[j] + tau2;
  }

  // Compute all log densities at once
  std::vector<double> logdens = dbivnorm_cpp(
      X, Y,
      mu_x, mu_y,
      sigma_xx, sigma_yy, sigma_xy,
      true  // log-scale
  );

  // Sum log-likelihood
  return std::accumulate(logdens.begin(), logdens.end(), 0.0);
}
