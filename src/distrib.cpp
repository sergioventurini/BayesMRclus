// distrib.cpp

#include "bayesmr.h"

// ------------------------------------------------------------
// Internal helpers (not visible outside this .cpp)
// ------------------------------------------------------------
namespace {

std::size_t recycle_size(std::size_t a,
                         std::size_t b,
                         std::size_t c = 1,
                         std::size_t d = 1){
  return std::max({a, b, c, d});
}

double recycle_get(const std::vector<double>& v, std::size_t i){
  return v[i % v.size()];
}

} // unnamed namespace

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

// Bivariate normal
std::vector<double> dbivnorm_cpp(const std::vector<double>& x_vec, const std::vector<double>& y_vec,
  const std::vector<double>& mu_x, const std::vector<double>& mu_y,
  const std::vector<double>& sigma_xx, const std::vector<double>& sigma_yy,
  const std::vector<double>& sigma_xy, bool logscale){
  const size_t n = x_vec.size();
  if (y_vec.size() != n || mu_x.size() != n || mu_y.size() != n ||
    sigma_xx.size() != n || sigma_yy.size() != n || sigma_xy.size() != n){
      Rcpp::stop("all input vectors must have the same length.");
  }

  std::vector<double> out(n);

  for (size_t i = 0; i < n; ++i){
    const double sxx = sigma_xx[i];
    const double syy = sigma_yy[i];
    const double sxy = sigma_xy[i];
    const double det = sxx * syy - sxy * sxy;

    // Invalid covariance?
    if (sxx <= 0.0 || syy <= 0.0 || det <= 0.0) {
      out[i] = logscale ?
        -std::numeric_limits<double>::infinity() :
        0.0;
      continue;
    }

    const double inv_det = 1.0 / det;
    const double dx = x_vec[i] - mu_x[i];
    const double dy = y_vec[i] - mu_y[i];

    const double q =
      inv_det * (syy * dx * dx - 2.0 * sxy * dx * dy + sxx * dy * dy);

    const double logdens =
      -0.5 * (log_2pi + std::log(det) + q);

    out[i] = logscale ? logdens : std::exp(logdens);
  }

  return out;
}

std::vector<double> dhalft(const std::vector<double>& x, const std::vector<double>& alpha,
  const std::vector<double>& nu, bool logscale){
  std::size_t N = recycle_size(x.size(), alpha.size(), nu.size());

  std::vector<double> out(N);

  for (double a : alpha)
    if (a <= 0.0)
      throw std::domain_error("The alpha (scale) parameter must be positive.");

  for (std::size_t i = 0; i < N; ++i) {
    double xi = recycle_get(x, i);
    double ai = recycle_get(alpha, i);
    double ni = recycle_get(nu, i);

    double log_dens =
        std::log(2.0)
        - std::log(ai)
        + R::lgammafn((ni + 1.0) / 2.0)
        - R::lgammafn(ni / 2.0)
        - 0.5 * std::log(M_PI * ni)
        - (ni + 1.0) / 2.0 * std::log(1.0 + (xi / ai) * (xi / ai) / ni);

    out[i] = logscale ? log_dens : std::exp(log_dens);
  }

  return out;
}

double dhalft_scalar(const double& x, const double& alpha,
  const double& nu, bool logscale){
  if (alpha <= 0.0)
    throw std::domain_error("The alpha (scale) parameter must be positive.");
  if (nu <= 0.0)
    throw std::domain_error("The nu (df) parameter must be positive.");

  // support
  if (x < 0.0)
    return logscale ? R_NegInf : 0.0;

  const double log_dens =
      std::log(2.0)
      - std::log(alpha)
      + R::lgammafn((nu + 1.0) / 2.0)
      - R::lgammafn(nu / 2.0)
      - 0.5 * std::log(M_PI * nu)
      - (nu + 1.0) / 2.0 *
        std::log(1.0 + (x / alpha) * (x / alpha) / nu);

  return logscale ? log_dens : std::exp(log_dens);
}

std::vector<double> pst(const std::vector<double>& q, const std::vector<double>& mu,
  const std::vector<double>& sigma, const std::vector<double>& nu,
  bool lower_tail, bool log_p){
  std::size_t N = recycle_size(q.size(), mu.size(), sigma.size(), nu.size());
  std::vector<double> out(N);

  for (double s : sigma)
    if (s <= 0.0)
      throw std::domain_error("The sigma parameter must be positive.");

  for (double n : nu)
    if (n <= 0.0)
      throw std::domain_error("The nu parameter must be positive.");

  for (std::size_t i = 0; i < N; ++i) {
    double qi = recycle_get(q, i);
    double mi = recycle_get(mu, i);
    double si = recycle_get(sigma, i);
    double ni = recycle_get(nu, i);

    if (ni > 1e6) {
      out[i] = R::pnorm(qi, mi, si, lower_tail, log_p);
    } else {
      double z = (qi - mi) / si;
      out[i] = R::pt(z, ni, lower_tail, log_p);
    }
  }

  return out;
}

std::vector<double> qst(const std::vector<double>& p, const std::vector<double>& mu,
  const std::vector<double>& sigma, const std::vector<double>& nu, bool lower_tail,
  bool log_p){
  std::size_t N = recycle_size(p.size(), mu.size(), sigma.size(), nu.size());
  std::vector<double> out(N);

  for (double s : sigma)
    if (s <= 0.0)
      throw std::domain_error("The sigma parameter must be positive.");

  for (double n : nu)
    if (n <= 0.0)
      throw std::domain_error("The nu parameter must be positive.");

  for (std::size_t i = 0; i < N; ++i) {
    double pi = recycle_get(p, i);
    double mi = recycle_get(mu, i);
    double si = recycle_get(sigma, i);
    double ni = recycle_get(nu, i);

    if (log_p)
      pi = std::exp(pi);

    if (pi < 0.0 || pi > 1.0)
      throw std::domain_error("p must be in [0,1].");

    if (ni > 1e6) {
      out[i] = R::qnorm(pi, mi, si, lower_tail, false);
    } else {
      out[i] = mi + si * R::qt(pi, ni, lower_tail, false);
    }
  }

  return out;
}

std::vector<double> qtrunc(const std::vector<double>& p, double a, double b,
  const std::vector<double>& mu, const std::vector<double>& sigma,
  const std::vector<double>& nu){
  if (a >= b)
    throw std::domain_error("Lower bound a is not less than upper bound b.");

  std::size_t N = recycle_size(p.size(), mu.size(), sigma.size(), nu.size());
  std::vector<double> out(N);

  for (std::size_t i = 0; i < N; ++i) {
    double pi = recycle_get(p, i);
    double mi = recycle_get(mu, i);
    double si = recycle_get(sigma, i);
    double ni = recycle_get(nu, i);

    if (pi < 0.0 || pi > 1.0)
      throw std::domain_error("p must be in [0,1].");

    double Ga = R::pt((a - mi) / si, ni, true, false);
    double Gb = R::pt((b - mi) / si, ni, true, false);

    out[i] = mi + si * R::qt(Ga + pi * (Gb - Ga), ni, true, false);
  }

  return out;
}

std::vector<double> rtrunc(std::size_t n, double a, double b,
  const std::vector<double>& mu, const std::vector<double>& sigma,
  const std::vector<double>& nu){
  if (a >= b)
    throw std::domain_error("Lower bound a is not less than upper bound b.");

  std::vector<double> u(n);
  for (std::size_t i = 0; i < n; ++i)
    u[i] = unif_rand();

  return qtrunc(u, a, b, mu, sigma, nu);
}

std::vector<double> rhalft(std::size_t n, const std::vector<double>& alpha, const std::vector<double>& nu){
  for (double a : alpha)
    if (a <= 0.0)
      throw std::domain_error("The alpha parameter must be positive.");

  std::vector<double> mu(1, 0.0);

  return rtrunc(n, 0.0, R_PosInf, mu, alpha, nu);
}
