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
  SEXP raiData,
  // SEXP raiD,
  // SEXP raix,
  // SEXP raing,
  SEXP radgamma,
  SEXP radbeta,
  // SEXP radalpha,
  SEXP rn,
  SEXP rp,
  SEXP rG,
  // SEXP rS,
  SEXP rtotiter,
  // SEXP radZ,
  SEXP rsigma2_beta,
  // SEXP rgamma_z,
  // SEXP reta,
  // SEXP rgamma_alpha,
  // SEXP rsigma2,
  // SEXP rlambda,
  SEXP rhyper_gammaj_gamma,
  SEXP rhyper_gammaj_psi2,
  SEXP rhyper_Gammaj_tau2,
  SEXP rhyper_gamma_mean,
  SEXP rhyper_gamma_var,
  SEXP rhyper_beta_mean,
  SEXP rhyper_beta_var,
  // SEXP rhyper_eta_a,
  // SEXP rhyper_eta_b,
  // SEXP rhyper_sigma2_a,
  // SEXP rhyper_sigma2_b,
  // SEXP rhyper_lambda,
  SEXP rverbose){
  SEXP rAns = NULL;
  int rAnsItems = 6;

  int n = Rf_asInteger(rn);
  int p = Rf_asInteger(rp);
  int G = Rf_asInteger(rG);
  // int S = Rf_asInteger(rS);
  int totiter = Rf_asInteger(rtotiter);
  // double gamma_z = Rf_asReal(rgamma_z);
  double sigma2_beta = Rf_asReal(rsigma2_beta);
  // double gamma_alpha = Rf_asReal(rgamma_alpha);
  double hyper_gammaj_gamma = Rf_asReal(rhyper_gammaj_gamma);
  double hyper_gammaj_psi2 = Rf_asReal(rhyper_gammaj_psi2);
  double hyper_Gammaj_tau2 = Rf_asReal(rhyper_Gammaj_tau2);
  double hyper_gamma_mean = Rf_asReal(rhyper_gamma_mean);
  double hyper_gamma_var = Rf_asReal(rhyper_gamma_var);
  double hyper_beta_mean = Rf_asReal(rhyper_beta_mean);
  double hyper_beta_var = Rf_asReal(rhyper_beta_var);
  // double hyper_sigma2_a = Rf_asReal(rhyper_sigma2_a);
  // double hyper_sigma2_b = Rf_asReal(rhyper_sigma2_b);
  int verbose = INTEGER(rverbose)[0];

  SEXP rgamma_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  SEXP rbeta_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*p*G)));
  // SEXP rz_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*n*p*G)));
  // SEXP ralpha_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  // SEXP reta_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  // SEXP rsigma2_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  // SEXP rlambda_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*G)));
  // SEXP rprob_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*S*G)));
  // SEXP rx_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*S)));
  // SEXP rx_ind_chain = PROTECT(Rf_allocVector(REALSXP, (totiter*S*G)));
  SEXP raccept = PROTECT(Rf_allocVector(REALSXP, G));
  SEXP rloglik = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogprior = PROTECT(Rf_allocVector(REALSXP, totiter));
  SEXP rlogpost = PROTECT(Rf_allocVector(REALSXP, totiter));

  bayesmr_mcmc_noclus(REAL(rgamma_chain), REAL(rbeta_chain), REAL(raccept), REAL(rloglik), REAL(rlogprior),
    REAL(rlogpost), REAL(raiData), REAL(radgamma), REAL(radbeta), REAL(rhyper_gammaj_gamma),
    REAL(rhyper_gammaj_psi2), REAL(rhyper_Gammaj_tau2), REAL(rhyper_gamma_mean), REAL(rhyper_gamma_var),
    REAL(rhyper_beta_mean), REAL(rhyper_beta_var), sigma2_beta, totiter, n, p, G, verbose);

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

// // //' Function for relabeling the parameter chain
// // //'
// // //' @param radtheta internal SEXP data structure
// // //' @param radz internal SEXP data structure
// // //' @param radeta internal SEXP data structure
// // //' @param radsigma2 internal SEXP data structure
// // //' @param radlambda internal SEXP data structure
// // //' @param radprob internal SEXP data structure
// // //' @param raix_ind internal SEXP data structure
// // //' @param rinit internal SEXP data structure
// // //' @param rM internal SEXP data structure
// // //' @param rR internal SEXP data structure
// // //'
// // //' @aliases bayesmr-internal
// // //' @aliases bayesmr_internal
// // //'
// // //' @rdname bayesmr_mcmc
// // //'
// // // [[Rcpp::export]] [[SV 20240116: commented to avoid errors in compilation]]
// RcppExport SEXP bayesmr_relabel(
//   SEXP radtheta,
//   SEXP radz,
//   SEXP radalpha,
//   SEXP radeta,
//   SEXP radsigma2,
//   SEXP radlambda,
//   SEXP radprob,
//   SEXP raix_ind,
//   SEXP rinit,
//   SEXP rn,
//   SEXP rp,
//   SEXP rS,
//   SEXP rM,
//   SEXP rR,
//   SEXP rG,
//   SEXP rverbose){
//   SEXP rAns = NULL;
//   const int rAnsItems = 8;

//   int init = Rf_asInteger(rinit);
//   int n = Rf_asInteger(rn);
//   int p = Rf_asInteger(rp);
//   int S = Rf_asInteger(rS);
//   int M = Rf_asInteger(rM);
//   int R = Rf_asInteger(rR);
//   int G = Rf_asInteger(rG);
//   int verbose = INTEGER(rverbose)[0];

//   relabel_celeux(REAL(radtheta), REAL(radz), REAL(radalpha), REAL(radeta), REAL(radsigma2), REAL(radlambda),
//       REAL(radprob), INTEGER(raix_ind), init, n, p, S, M, R, G, verbose);

//   // packing results
//   PROTECT(rAns = Rf_allocVector(VECSXP, rAnsItems));

//   SET_VECTOR_ELT(rAns, 0, radtheta);
//   SET_VECTOR_ELT(rAns, 1, radz);
//   SET_VECTOR_ELT(rAns, 2, radalpha);
//   SET_VECTOR_ELT(rAns, 3, radeta);
//   SET_VECTOR_ELT(rAns, 4, radsigma2);
//   SET_VECTOR_ELT(rAns, 5, radlambda);
//   SET_VECTOR_ELT(rAns, 6, radprob);
//   SET_VECTOR_ELT(rAns, 7, raix_ind);

//   // cleanup and return
//   UNPROTECT(1);  // rAns

//   return rAns;
// }

// // //' Packing of the parameter chain to be run before relabeling
// // //'
// // //' @rdname bayesmr_mcmc
// // //'
// // // [[Rcpp::export]] [[SV 20240116: commented to avoid errors in compilation]]
// RcppExport SEXP bayesmr_pack_par(
//   SEXP radz,
//   SEXP radalpha,
//   SEXP radlambda,
//   SEXP rn,
//   SEXP rp,
//   SEXP rM,
//   SEXP rG){
//   int n = Rf_asInteger(rn);
//   int p = Rf_asInteger(rp);
//   int M = Rf_asInteger(rM);
//   int G = Rf_asInteger(rG);
//   int r = n*(n - 1)/2;

//   SEXP rAns = PROTECT(Rf_allocVector(REALSXP, M*(r + 1)*G));

//   pack_par(REAL(rAns), REAL(radz), REAL(radalpha), REAL(radlambda), n, p, M, G);

//   // cleanup and return
//   UNPROTECT(1);  // rAns

//   return rAns;
// }
