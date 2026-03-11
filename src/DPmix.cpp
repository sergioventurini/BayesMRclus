// DPmix.cpp

#include "bayesmr.h"

// Full conditional of gamma
double full_cond_gamma(
  const double mu_gamma,
  const double var_gamma,
  const std::vector<double>& gamma_hat,
  const std::vector<double>& Gamma_hat,
  const std::vector<double>& sigma2_X,
  const std::vector<double>& sigma2_Y,
  const double psi,
  const double tau,
  const std::vector<double>& beta,
  const std::vector<int>& xi)
{
  const int n = gamma_hat.size();
  const double var_inv = 1.0 / var_gamma;

  double A = var_inv;
  double B = mu_gamma * var_inv;

  for (int j = 0; j < n; ++j) {
    const double b = beta[xi[j]];
    const double s11 = psi*psi + sigma2_X[j];
    const double s12 = b * psi*psi;
    const double s21 = s12;
    const double s22 = b*b*psi*psi + tau*tau + sigma2_Y[j];

    double inv11, inv12, inv21, inv22;
    inverse_2x2(s11, s12, s21, s22, inv11, inv12, inv21, inv22);

    // A_j' Σ^-1 A_j with A_j = [1, b]
    const double A_term = inv11 + b*(inv12 + inv21) + b*b*inv22;
    A += A_term;

    // A_j' Σ^-1 y_j,   y_j = [γ_hat_j, Γ_hat_j]
    const double B_term_1 = inv11*gamma_hat[j] + inv12*Gamma_hat[j];
    const double B_term_2 = inv21*gamma_hat[j] + inv22*Gamma_hat[j];
    B += B_term_1 + b*B_term_2;
  }

  const double mean = B / A;
  const double sd = std::sqrt(1.0 / A);

  return R::rnorm(mean, sd);
}

// Full conditional of beta (log density)
double full_cond_beta(double beta, double mu_beta, double sigma2_beta,
  const std::vector<double>& gamma_hat,
  const std::vector<double>& Gamma_hat,
  const std::vector<double>& sigma2X,
  const std::vector<double>& sigma2Y,
  double gamma, double psi, double tau) {

  double log_prior = -0.5 / sigma2_beta * (beta - mu_beta) * (beta - mu_beta);

  int n = gamma_hat.size();
  double log_lik = bayesmr_logLik(
    beta, gamma, psi*psi, tau*tau, n,
    gamma_hat.data(), Gamma_hat.data(),
    sigma2X.data(), sigma2Y.data()
  );

  return log_prior + log_lik;
}

// Joint full conditional of psi and tau (log density)
double full_cond_psitau(const std::vector<double>& beta,
  double alpha_psi, double nu_psi,
  double alpha_tau, double nu_tau,
  const std::vector<double>& gamma_hat,
  const std::vector<double>& Gamma_hat,
  const std::vector<double>& sigma2X,
  const std::vector<double>& sigma2Y,
  double gamma, double psi, double tau) {

  double log_prior_psi = dhalft_scalar(psi, alpha_psi, nu_psi, true);
  double log_prior_tau = dhalft_scalar(tau, alpha_tau, nu_tau, true);

  int n = gamma_hat.size();
  double log_lik = bayesmr_logLik_vec(
    beta, gamma, psi*psi, tau*tau, n,
    gamma_hat.data(), Gamma_hat.data(),
    sigma2X.data(), sigma2Y.data()
  );

  return log_lik + log_prior_psi + log_prior_tau;
}

// Full conditional of log(psi) (log density)
double full_cond_psi(std::vector<double> beta, double alpha_psi, double nu_psi,
  const std::vector<double>& gamma_hat,
  const std::vector<double>& Gamma_hat,
  const std::vector<double>& sigma2X,
  const std::vector<double>& sigma2Y,
  double gamma, double psi, double tau) {

  double log_prior = dhalft_scalar(psi, alpha_psi, nu_psi, true) + std::log(psi);

  int n = gamma_hat.size();
  double log_lik = bayesmr_logLik_vec(
    beta, gamma, psi*psi, tau*tau, n,
    gamma_hat.data(), Gamma_hat.data(),
    sigma2X.data(), sigma2Y.data()
  );

  return log_prior + log_lik;
}

// Full conditional of log(tau) (log density)
double full_cond_tau(std::vector<double> beta, double alpha_tau, double nu_tau,
  const std::vector<double>& gamma_hat,
  const std::vector<double>& Gamma_hat,
  const std::vector<double>& sigma2X,
  const std::vector<double>& sigma2Y,
  double gamma, double psi, double tau) {

  double log_prior = dhalft_scalar(tau, alpha_tau, nu_tau, true) + std::log(tau);

  int n = gamma_hat.size();
  double log_lik = bayesmr_logLik_vec(
    beta, gamma, psi*psi, tau*tau, n,
    gamma_hat.data(), Gamma_hat.data(),
    sigma2X.data(), sigma2Y.data()
  );

  return log_prior + log_lik;
}

// Normalize weights
std::vector<double> normalize_weights(const std::vector<double>& prob) {
  std::vector<double> weight(prob.size());

  // Find max of non-inf values
  double const_val = -std::numeric_limits<double>::infinity();
  for (double p : prob) {
    if (std::isfinite(p) && p > const_val) {
      const_val = p;
    }
  }

  // Compute exp and sum
  double exp_sum = 0.0;
  for (size_t i = 0; i < prob.size(); ++i) {
    if (std::isfinite(prob[i])) {
      weight[i] = std::exp(prob[i] - const_val);
      exp_sum += weight[i];
    } else {
      weight[i] = 0.0;
    }
  }

  // Normalize
  if (exp_sum > 0.0) {
    for (double& w : weight) {
      w /= exp_sum;
    }
  }

  return weight;
}

// Count occurrences excluding index j
int count_excluding(const std::vector<int>& xi, int k, size_t exclude_idx) {
  int count = 0;
  for (size_t i = 0; i < xi.size(); ++i) {
    if (i != exclude_idx && xi[i] == k) {
      ++count;
    }
  }
  return count;
}

// Relabel xi to be contiguous starting from 0
void relabel_xi(std::vector<int>& xi) {
  std::set<int> unique_vals(xi.begin(), xi.end());
  std::map<int, int> old_to_new;
  int new_label = 0;
  for (int val : unique_vals) {
    old_to_new[val] = new_label++;
  }
  for (int& val : xi) {
    val = old_to_new[val];
  }
}
