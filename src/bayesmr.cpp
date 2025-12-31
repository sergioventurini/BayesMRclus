// bayesmr.cpp

#include "bayesmr.h"

// Note: RcppExport is an alias for extern "C"

// [[SV 20240702: documentation and exporting has been commented to avoid
//                errors in compilation]]

// //' Internal functions for MCMC simulation.
// //'
// //' For internal use only.
// //'
// //' @param raiD internal SEXP data structure
// //' @param raix internal SEXP data structure
// //' @param raing internal SEXP data structure
// //' @param radalpha internal SEXP data structure
// //' @param rn internal SEXP data structure
// //' @param rp internal SEXP data structure
// //' @param rG internal SEXP data structure
// //' @param rS internal SEXP data structure
// //' @param rtotiter internal SEXP data structure
// //' @param radZ internal SEXP data structure
// //' @param rgamma_z internal SEXP data structure
// //' @param reta internal SEXP data structure
// //' @param rgamma_alpha internal SEXP data structure
// //' @param rsigma2 internal SEXP data structure
// //' @param rlambda internal SEXP data structure
// //' @param rhyper_eta_a internal SEXP data structure
// //' @param rhyper_eta_b internal SEXP data structure
// //' @param rhyper_sigma2_a internal SEXP data structure
// //' @param rhyper_sigma2_b internal SEXP data structure
// //' @param rhyper_lambda internal SEXP data structure
// //' @param rfamily internal SEXP data structure
// //' @param rverbose internal SEXP data structure
// //'
// //' @aliases bayesmr-internal
// //' @aliases bayesmr_internal
// //'
// // [[Rcpp::export]]
RcppExport SEXP bayesmr_mcmc(
  SEXP radData,
  SEXP radgamma,
  SEXP radbeta,
  SEXP rn,
  SEXP rp,
  SEXP rG,
  SEXP rtotiter,
  SEXP rsigma2_beta,
  SEXP rhyper_gammaj_psi2,
  SEXP rhyper_Gammaj_tau2,
  SEXP rhyper_gamma_mean,
  SEXP rhyper_gamma_var,
  SEXP rhyper_beta_mean,
  SEXP rhyper_beta_var,
  SEXP rverbose){
  SEXP rAns = NULL;
  int rAnsItems = 6;

  int n = Rf_asInteger(rn);
  int p = Rf_asInteger(rp);
  int G = Rf_asInteger(rG);
  int totiter = Rf_asInteger(rtotiter);
  double sigma2_beta = Rf_asReal(rsigma2_beta);
  double gamma_val = Rf_asReal(radgamma);
  double beta_val = Rf_asReal(radbeta);
  double hyper_gammaj_psi2 = Rf_asReal(rhyper_gammaj_psi2);
  double hyper_Gammaj_tau2 = Rf_asReal(rhyper_Gammaj_tau2);
  double hyper_gamma_mean = Rf_asReal(rhyper_gamma_mean);
  double hyper_gamma_var = Rf_asReal(rhyper_gamma_var);
  double hyper_beta_mean = Rf_asReal(rhyper_beta_mean);
  double hyper_beta_var = Rf_asReal(rhyper_beta_var);
  int verbose = INTEGER(rverbose)[0];

  SEXP rgamma_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP rbeta_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP raccept = PROTECT(Rf_allocVector(REALSXP, G));
  SEXP rloglik = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogprior = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogpost = PROTECT(Rf_allocVector(REALSXP, totiter));

  bayesmr_mcmc_noclus(REAL(rgamma_chain), REAL(rbeta_chain), REAL(raccept),
    REAL(rloglik), REAL(rlogprior), REAL(rlogpost), REAL(radData), gamma_val,
    beta_val, hyper_gammaj_psi2, hyper_Gammaj_tau2, hyper_gamma_mean, hyper_gamma_var,
    hyper_beta_mean, hyper_beta_var, sigma2_beta, totiter, n, p, G, verbose);

  // packing results
  PROTECT(rAns = Rf_allocVector(VECSXP, rAnsItems));

  SET_VECTOR_ELT(rAns, 0, rgamma_chain);
  SET_VECTOR_ELT(rAns, 1, rbeta_chain);
  SET_VECTOR_ELT(rAns, 2, raccept);
  SET_VECTOR_ELT(rAns, 3, rloglik);
  SET_VECTOR_ELT(rAns, 4, rlogprior);
  SET_VECTOR_ELT(rAns, 5, rlogpost);

  // cleanup and return
  UNPROTECT(7);  // rAns, rgamma_chain, rbeta_chain, raccept, rloglik, rlogprior, rlogpost

  return rAns;
}

// //' Internal functions for MCMC simulation.
// //'
// //' For internal use only.
// //'
// //' @param raiD internal SEXP data structure
// //' @param raix internal SEXP data structure
// //' @param raing internal SEXP data structure
// //' @param radalpha internal SEXP data structure
// //' @param rn internal SEXP data structure
// //' @param rp internal SEXP data structure
// //' @param rG internal SEXP data structure
// //' @param rS internal SEXP data structure
// //' @param rtotiter internal SEXP data structure
// //' @param radZ internal SEXP data structure
// //' @param rgamma_z internal SEXP data structure
// //' @param reta internal SEXP data structure
// //' @param rgamma_alpha internal SEXP data structure
// //' @param rsigma2 internal SEXP data structure
// //' @param rlambda internal SEXP data structure
// //' @param rhyper_eta_a internal SEXP data structure
// //' @param rhyper_eta_b internal SEXP data structure
// //' @param rhyper_sigma2_a internal SEXP data structure
// //' @param rhyper_sigma2_b internal SEXP data structure
// //' @param rhyper_lambda internal SEXP data structure
// //' @param rfamily internal SEXP data structure
// //' @param rverbose internal SEXP data structure
// //'
// //' @aliases bayesmr-internal
// //' @aliases bayesmr_internal
// //'
// // [[Rcpp::export]]
RcppExport SEXP bayesmr_mcmc_het(
  SEXP radData,
  SEXP radgamma,
  SEXP radbeta,
  SEXP radpsi,
  SEXP radtau,
  SEXP rn,
  SEXP rp,
  SEXP rG,
  SEXP rtotiter,
  SEXP rsigma2_beta,
  SEXP rC_psi,
  SEXP rC_tau,
  SEXP rhyper_gamma_mean,
  SEXP rhyper_gamma_var,
  SEXP rhyper_beta_mean,
  SEXP rhyper_beta_var,
  SEXP rhyper_psi_alpha,
  SEXP rhyper_psi_nu,
  SEXP rhyper_tau_alpha,
  SEXP rhyper_tau_nu,
  SEXP rverbose){
  SEXP rAns = NULL;
  int rAnsItems = 8;

  int n = Rf_asInteger(rn);
  int p = Rf_asInteger(rp);
  int G = Rf_asInteger(rG);
  int totiter = Rf_asInteger(rtotiter);
  double sigma2_beta = Rf_asReal(rsigma2_beta);
  double C_psi = Rf_asReal(rC_psi);
  double C_tau = Rf_asReal(rC_tau);
  double gamma_val = Rf_asReal(radgamma);
  double beta_val = Rf_asReal(radbeta);
  double psi_val = Rf_asReal(radpsi);
  double tau_val = Rf_asReal(radtau);
  double hyper_gamma_mean = Rf_asReal(rhyper_gamma_mean);
  double hyper_gamma_var = Rf_asReal(rhyper_gamma_var);
  double hyper_beta_mean = Rf_asReal(rhyper_beta_mean);
  double hyper_beta_var = Rf_asReal(rhyper_beta_var);
  double hyper_psi_alpha = Rf_asReal(rhyper_psi_alpha);
  double hyper_psi_nu = Rf_asReal(rhyper_psi_nu);
  double hyper_tau_alpha = Rf_asReal(rhyper_tau_alpha);
  double hyper_tau_nu = Rf_asReal(rhyper_tau_nu);
  int verbose = INTEGER(rverbose)[0];

  SEXP rgamma_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP rbeta_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP rpsi_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP rtau_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP raccept = PROTECT(Rf_allocVector(REALSXP, 3*G));
  SEXP rloglik = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogprior = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogpost = PROTECT(Rf_allocVector(REALSXP, totiter));

  bayesmr_mcmc_noclus_het(REAL(rgamma_chain), REAL(rbeta_chain),
    REAL(rpsi_chain), REAL(rtau_chain), REAL(raccept),
    REAL(rloglik), REAL(rlogprior), REAL(rlogpost), REAL(radData), gamma_val,
    beta_val, psi_val, tau_val, hyper_psi_alpha, hyper_psi_nu, hyper_tau_alpha,
    hyper_tau_nu, hyper_gamma_mean, hyper_gamma_var, hyper_beta_mean, hyper_beta_var,
    sigma2_beta, C_psi, C_tau, totiter, n, p, G, verbose);

  // packing results
  PROTECT(rAns = Rf_allocVector(VECSXP, rAnsItems));

  SET_VECTOR_ELT(rAns, 0, rgamma_chain);
  SET_VECTOR_ELT(rAns, 1, rbeta_chain);
  SET_VECTOR_ELT(rAns, 2, rpsi_chain);
  SET_VECTOR_ELT(rAns, 3, rtau_chain);
  SET_VECTOR_ELT(rAns, 4, raccept);
  SET_VECTOR_ELT(rAns, 5, rloglik);
  SET_VECTOR_ELT(rAns, 6, rlogprior);
  SET_VECTOR_ELT(rAns, 7, rlogpost);

  // cleanup and return
  UNPROTECT(9);  // rAns, rgamma_chain, rbeta_chain, rpsi_chain, rtau_chain, raccept, rloglik, rlogprior, rlogpost

  return rAns;
}
