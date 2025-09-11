// bayesmr_post.cpp

#include "bayesmr.h"

// full conditional distribution of the beta parameter
void logpost_beta(double* lpost, const double beta, const double gamma, const double mu_beta,
  const double sigma2_beta, const double psi2, const double tau2, int n,
  const double* Gammahat_j, const double* sigma2_Y){
  double h2_j = 0, tmp = 0;
  for(int j = 0; j < n; j++){
    h2_j = pow(beta, 2)*psi2 + sigma2_Y[j] + tau2;
    tmp += log(h2_j) + pow(Gammahat_j[j] - beta*gamma, 2)/h2_j;
  }

  *lpost = -0.5*(tmp + pow(beta - mu_beta, 2)/sigma2_beta);
}
