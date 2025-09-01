// bayesmr_post.cpp

#include "bayesmr.h"

// full conditional distribution of the beta parameter
void logpost_beta(double* lpost, const double beta, const double gamma, const double mu_beta,
  const double sigma2_beta, const double psi2, const double tau2, int n,
  const double* Gammahat_j, const double* sigma2_Y){
  double* tau2_j = new double[n];
  double* h2_j = new double[n];
  double tmp = 0;
  for(int i = 0; i < n; i++){
    tau2_j[i] = sigma2_Y[i] + tau2;
    h2_j[i] = pow(beta, 2)*psi2 + tau2_j[i];
    tmp += log(h2_j[i]) + pow(Gammahat_j[i] - beta*gamma, 2)/h2_j[i];
  }
  double loglik_Gammahat = -0.5*tmp;
  double logprior_beta = -0.5*pow(beta - mu_beta, 2)/sigma2_beta;

  *lpost = loglik_Gammahat + logprior_beta;

  delete[] tau2_j;
  delete[] h2_j;
}
