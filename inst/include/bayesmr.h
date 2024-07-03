#ifndef BAYESMR_H
#define BAYESMR_H

// #define ARMA_NO_DEBUG // this macro disables bounds checks in the Armadillo
                         // library making the code faster (but also more
                         // frail!); the suggestion is to disable it only for
                         // the final release (see
                         // http://arma.sourceforge.net/docs.html)

#include <R.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <climits>
#include <string>
#include <algorithm>

static const double machine_eps = 2.220446049250313080847e-16;
static const double log_pi = std::log(M_PI);
static const double log_2pi = std::log(2.0 * M_PI);
static const double log_two = std::log(2.0);

// MAIN FUNCTIONS -----------------------------------------------------------------------------------------------------
RcppExport SEXP bayesmr_mcmc(SEXP radData, SEXP radgamma, SEXP radbeta,
  SEXP rn, SEXP rp, SEXP rG, SEXP rtotiter, SEXP rsigma2_beta,
  SEXP rhyper_gammaj_gamma, SEXP rhyper_gammaj_psi2, SEXP rhyper_Gammaj_tau2,
  SEXP rhyper_gamma_mean, SEXP rhyper_gamma_var, SEXP rhyper_beta_mean,
  SEXP rhyper_beta_var, SEXP rverbose);
// RcppExport SEXP bayesmr_relabel(SEXP radtheta, SEXP radz, SEXP radalpha,
//   SEXP radeta, SEXP radsigma2, SEXP radlambda, SEXP radprob, SEXP raix_ind,
//   SEXP rinit, SEXP rn, SEXP rp, SEXP rS, SEXP rM, SEXP rR, SEXP rG,
//   SEXP rverbose);
// RcppExport SEXP bayesmr_pack_par(SEXP radz, SEXP radalpha, SEXP radlambda,
//   SEXP rn, SEXP rp, SEXP rM, SEXP rG);

// DISTRIBUTION FUNCTIONS ---------------------------------------------------------------------------------------------
void dprodber(double* prob, const int* d, const double* pi, int m,
  int logscale);
void dmultinorm(double* dens, const double* x, const double* mean,
  const double* sigma, int n, int p, int logscale = 0);
void rmultinorm(double* dev, int n, const double* mean,
  const double* sigma, int p);
void dinvgamma(double* dens, const double* x, const double alpha,
  const double beta, int n, int logscale);
void rinvgamma(double* dev, int n, const double alpha, const double beta);
void ddirichlet(double* dens, const double* x, const double* par, int n, int p,
  int logscale);
void rdirichlet(double* dev, int n, const double* par, int p);
void loglik_rbmds_binom(double* loglik, const int* d, const double* z,
  const double alpha, int n, int p, int S);
void loglik_bayesmr(double* loglik, const int* d, const double* z,
  const double* alpha, const double* sigma2, const double* lambda,
  const int* x, int n, int p, int S, int G, const char* family);

// MATRIX UTILITIES ---------------------------------------------------------------------------------------------------
void colsums(double* colsums, const double* A, int nrows, int ncols);
void rowsums(double* rowsums, const double* A, int nrows, int ncols);
arma::vec dissM2V(const arma::mat& d);
arma::vec mat2vec(const arma::mat& A, const int& j);
arma::mat vec2mat(const arma::mat& A, const int& j, const arma::vec& v);
arma::vec mahalanobis(const arma::mat& x, const arma::vec& center,
  const arma::mat& cov);

// MCMC SIMULATION ----------------------------------------------------------------------------------------------------
void bayesmr_mcmc_noclus(double* gamma_chain, double* beta_chain, double* accept,
  double* loglik, double* logprior, double* logpost, double* data, double* gamma,
  double* beta, const double* hyper_gammaj_gamma, const double* hyper_gammaj_psi2,
  const double* rhyper_Gammaj_tau2, const double* rhyper_gamma_mean,
  const double* rhyper_gamma_var, const double* rhyper_beta_mean,
  const double* rhyper_beta_var, const double sigma2_beta,
  int totiter, int n, int p, int G, int verbose);

// RELABEL ALGORITHMS -------------------------------------------------------------------------------------------------
// void relabel_celeux(double* theta, double* z_chain, double* alpha_chain,
//   double* eta_chain, double* sigma2_chain, double* lambda_chain,
//   double* prob_chain, int* x_ind_chain, int init, int n, int p, int S, int M,
//   int R, int G, int verbose);
// void pack_par(double* theta, const double* z, const double* alpha,
//   const double* lambda, int n, int p, int M, int G);

// UTILITIES ----------------------------------------------------------------------------------------------------------
void logit(double* res, const double* p, int n);
void expit(double* res, const double* x, int n);
void exp_vec(double* res, const double* x, int n);
double euclidean(const double *x, int nr, int nc, int i1, int i2);
void dist(double* d, const double* x, int nr, int nc);
bool any_na_nan(const arma::vec x, const int& n);
void sample_no_rep(int n, double* p, int* perm, int nans, int* ans);
void tableC(int* counts, const int* x, int nelem, int ndistelem);
int factorial(const int& x);
void permutations(int* perm, int n, int nperm, int byrow);
void which_min(int* ans, const double* r, int n);

// For registration
void R_init_bayesmr(DllInfo *dll);

#endif
