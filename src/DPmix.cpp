// DPmix.cpp

#include "bayesmr.h"

// Full conditional of gamma
double full_cond_gamma(
  const double mu_gamma,
  const double var_gamma,
  const std::vector<double>& gamma_hat,
  const std::vector<double>& Gamma_hat,
  const std::vector<double>& psi2_j,    // sigma2_X + rhyper_gammaj_psi2
  const std::vector<double>& tau2_j,    // sigma2_Y + rhyper_Gammaj_tau2
  const std::vector<double>& beta,
  const std::vector<int>& xi)
{
  const int n = gamma_hat.size();
  const double var_inv = 1.0 / var_gamma;

  double A = var_inv;
  double B = mu_gamma * var_inv;

  for (int j = 0; j < n; ++j) {
    const double b = beta[xi[j]];
    const double s11 = psi2_j[j];
    const double s12 = b * psi2_j[j];
    const double s21 = s12;
    const double s22 = b*b*psi2_j[j] + tau2_j[j];

    double inv11, inv12, inv21, inv22;
    inverse_2x2(s11, s12, s21, s22,
      inv11, inv12, inv21, inv22);

    // A_j' Σ^-1 A_j with A_j = [1, b]
    const double AjSinvAj = inv11 + b*(inv12 + inv21) + b*b*inv22;
    A += AjSinvAj;

    // A_j' Σ^-1 y_j,   y_j = [γ_hat_j, Γ_hat_j]
    const double Sinv_y1 = inv11*gamma_hat[j] + inv12*Gamma_hat[j];
    const double Sinv_y2 = inv21*gamma_hat[j] + inv22*Gamma_hat[j];
    B += Sinv_y1 + b*Sinv_y2;
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
  double gamma, double psi2) {

  double log_prior = -0.5 / sigma2_beta * (beta - mu_beta) * (beta - mu_beta);

  int n = gamma_hat.size();
  double log_lik = bayesmr_logLik(
    beta, gamma, psi2, 0.0, n,
    gamma_hat.data(), Gamma_hat.data(),
    sigma2X.data(), sigma2Y.data()
  );

  return log_prior + log_lik;
}

// Normalize weights
std::vector<double> normalize_weights(const std::vector<double>& prob) {
  std::vector<double> weight(prob.size());

  // Find mean of non-inf values
  double sum = 0.0;
  int count = 0;
  for (double p : prob) {
    if (std::isfinite(p)) {
      sum += p;
      ++count;
    }
  }
  double const_val = (count > 0) ? sum / count : 0.0;

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
