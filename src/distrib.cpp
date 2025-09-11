// distrib.cpp

#include "bayesmr.h"

// Probability function for the product of independent bernoulli random variables
void dprodber(double* prob, const int* d, const double* pi, int m, bool logscale = true){
  double prob_tmp = 0;
  if(!logscale){
    *prob = 0;
    for(int i = 0; i < m; i++){
      prob_tmp = d[i]*log(pi[i]) + (1 - d[i])*log(1 - pi[i]);
      if(ISNAN(prob_tmp)){
        prob_tmp = log(pow(pi[i], d[i])*pow(1 - pi[i], 1 - d[i]));
      }
      *prob += prob_tmp;
    }
  }
  else{
    *prob = 1;  
    for(int i = 0; i < m; i++){
      prob_tmp = pow(pi[i], d[i])*pow(1 - pi[i], 1 - d[i]);
      *prob *= prob_tmp;
    }
  }
}

// Density function for a multivariate normal random variable
arma::vec dmvnorm_fast(const arma::mat& X,
                       const arma::vec& mu,
                       const arma::mat& Sigma,
                       bool logscale = true){
  
  int n = X.n_rows;
  int k = X.n_cols;
  
  // Pre-compute Cholesky and constants
  arma::mat L = arma::chol(Sigma, "lower");
  double log_det_sigma = 2.0 * arma::sum(arma::log(L.diag()));
  
  arma::vec result(n);
  
  for(int i = 0; i < n; i++){
    arma::vec diff = X.row(i).t() - mu;
    arma::vec z = arma::solve(arma::trimatl(L), diff);
    
    result(i) = -0.5 * (k * log_2pi + log_det_sigma + arma::dot(z, z));
  }
  
  if(!logscale){
    result = arma::exp(result);
  }
  
  return result;
}

// Density function for a multivariate normal random variable
void dmultinorm(double* log_dens, const double* x, const double* mean, const double* sigma, int n, int p,
                bool logscale = true){
  arma::mat x_arma(n, p);
  arma::vec mean_arma(p);
  arma::mat sigma_arma(p, p);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      x_arma(i, j) = x[n*j + i];
    }
    mean_arma(j) = mean[j];
    for(int k = 0; k < p; k++){
      sigma_arma(j, k) = sigma[p*k + j];
    }
  }

  double logdet, sign;
  arma::log_det(logdet, sign, sigma_arma);
  arma::vec distval = mahalanobis(x_arma, mean_arma, sigma_arma);
  for(int i = 0; i < n; i++){
    log_dens[i] = -0.5*(p * log(2.0*M_PI) + logdet + distval[i]);
  }

  if(!logscale){
    for(int i = 0; i < n; i++){
      log_dens[i] = exp(log_dens[i]);
    }
  }
}

// Multivariate normal random deviates
void rmultinorm(double* dev, int n, const double* mean, const double* sigma, int p){
  arma::vec mean_arma(p);
  arma::mat sigma_arma(p, p);
  for(int j = 0; j < p; j++){
    mean_arma(j) = mean[j];
    for(int i = 0; i < p; i++){
      sigma_arma(i, j) = sigma[i + p*j];
    }
  }

  arma::mat V(p, p);
  arma::mat U(p, p);
  arma::vec d(p);
  arma::svd(U, d, V, sigma_arma);
  arma::mat U_trans = trans(U);

  arma::mat tmp(p, p);
  for(int i = 0; i < p; i++){
    for(int j = 0; j < p; j++){
      tmp(i, j) = U_trans(i, j)*sqrt(d(i));
    }
  }
  arma::mat sigmaSVDsqrt = trans(V * tmp);
  arma::mat normdev(n, p);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      normdev(i, j) = R::rnorm(0, 1);
    }
  }
  arma::mat out_tmp = normdev * sigmaSVDsqrt;
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      dev[n*j + i] = out_tmp(i, j) + mean_arma(j);
    }
  }
}

// Density function for an inverse gamma random variable
void dinvgamma(double* log_dens, const double* x, const double alpha, const double beta, int n,
               bool logscale = true){
  if((alpha <= 0) || (beta <= 0)){
    error("alpha (shape) and beta (scale) parameters in dinvgamma() need to be both strictly positive.\n");
  }

  double lbeta = log(beta);
  double lgalpha = R::lgammafn(alpha);
  for(int i = 0; i < n; i++){
    log_dens[i] = alpha * lbeta - lgalpha - (alpha + 1) * log(x[i]) - (beta/x[i]);
  }

  if(!logscale){
    for(int i = 0; i < n; i++){
      log_dens[i] = exp(log_dens[i]);
    }
  }
}

// Inverse gamma random deviates
void rinvgamma(double* dev, int n, const double alpha, const double beta){
  if((alpha <= 0) || (beta <= 0)){
    error("alpha (shape) and beta (scale) parameters in rinvgamma() need to be both strictly positive.\n");
  }

  for(int i = 0; i < n; i++){
    dev[i] = 1.0/R::rgamma(alpha, 1.0/beta);
  }
}

// Density function for a Dirichlet random variable
void ddirichlet(double* log_dens, const double* x, const double* par, int n, int p,
                bool logscale = true){
  double tmp = 0;
  arma::mat x_arma(n, p);
  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      x_arma(i, j) = x[n*j + i];
    }
  }

  tmp = 0;
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      tmp += x_arma(i, j);
      if(x_arma(i, j) < 0 || x_arma(i, j) > 1){
        error("some elements of x outside the [0, 1] range.\n");
      }
    }
    if(std::fabs(tmp - 1.0) >= sqrt(std::numeric_limits<double>::epsilon())){
      error("some rows of x sum to a value different from 1.\n");
    }
  }

  double logD = 0, salpha = 0, s = 0;
  for(int j = 0; j < p; j++){
    logD += R::lgammafn(par[j]);
    salpha += par[j];
  }
  logD = logD - R::lgammafn(salpha);
  for(int i = 0; i < n; i++){
    s = 0;
    for(int j = 0; j < p; j++){
      s += (par[j] - 1)*log(x_arma(i, j));
    }
    log_dens[i] = s - logD;
    if(!logscale){
      log_dens[i] = exp(log_dens[i]);
    }
  }
}

// Dirichlet random deviates
void rdirichlet(double* dev, int n, const double* par, int p){
  arma::mat z(n, p, arma::fill::zeros);
  arma::vec s(n, arma::fill::zeros);

  for(int j = 0; j < p; j++){
    for(int i = 0; i < n; i++){
      z(i, j) = R::rgamma(par[j], 1.0);
      s(i) += z(i, j);
    }
  }
  for(int i = 0; i < n; i++){
    for(int j = 0; j < p; j++){
      dev[n*j + i] = z(i, j)/s(i);
    }
  }
}

bool bivnorm_validate_covariance(double sigma_xx, double sigma_yy, double sigma_xy){
    if(sigma_xx <= 0.0 || sigma_yy <= 0.0) return false;
    double det = sigma_xx * sigma_yy - sigma_xy * sigma_xy;
    return det > 0.0;
}

// Bivariate normal
std::vector<double> dbivnorm_cpp(const std::vector<double>& x_vec,
                             const std::vector<double>& y_vec,
                             double mu_x, double mu_y,
                             double sigma_xx, double sigma_yy, double sigma_xy,
                             bool logscale){
  size_t n = x_vec.size();
  if(n != y_vec.size()){
    // Handle size mismatch - return empty vector or throw exception
    return std::vector<double>();
  }
  
  std::vector<double> result(n);
  
  // Validate covariance matrix
  if(!bivnorm_validate_covariance(sigma_xx, sigma_yy, sigma_xy)){
    std::fill(result.begin(), result.end(), logscale ? neg_inf : 0.0);
    return result;
  }
  
  // Precompute constants
  double det_sigma = sigma_xx * sigma_yy - sigma_xy * sigma_xy;
  double inv_det = 1.0 / det_sigma;
  double log_const = -0.5 * (log_2pi + std::log(det_sigma));
  
  // Compute density for each observation
  for(size_t i = 0; i < n; ++i){
    double dx = x_vec[i] - mu_x;
    double dy = y_vec[i] - mu_y;
    
    double quad_form = inv_det * (sigma_yy * dx * dx - 2.0 * sigma_xy * dx * dy + sigma_xx * dy * dy);
    double log_dens = log_const - 0.5 * quad_form;
    
    result[i] = logscale ? log_dens : std::exp(log_dens);
  }
  
  return result;
}
