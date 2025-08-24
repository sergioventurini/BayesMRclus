// mcmc.cpp

#include "bayesmr.h"

void bayesmr_mcmc_noclus(
  double* gamma_chain,
  double* beta_chain,
  double* accept,
  double* loglik,
  double* logprior,
  double* logpost,
  double* data,
  double gamma_p,                             // gamma starting value
  double beta_p,                              // beta starting value
  const double hyper_gammaj_psi2,
  const double rhyper_Gammaj_tau2,
  const double rhyper_gamma_mean,
  const double rhyper_gamma_var,
  const double rhyper_beta_mean,
  const double rhyper_beta_var,
  const double sigma2_beta,                   // beta proposal spread
  int totiter,
  int n,
  int p,
  int G,
  int verbose){
  int niter = 0;

  double beta_prop = 0, beta_old = 0, accept_beta = 0;

  double* gammahat_j = new double[n];
  double* Gammahat_j = new double[n];
  double* sigma2_X = new double[n];
  double* sigma2_Y = new double[n];
  double* psi2_j = new double[n];
  double* tau2_j = new double[n];
  double* h2_j = new double[n];
  double sum_A_beta_1 = 0, sum_A_beta_2 = 0, sum_B_beta_1 = 0, sum_B_beta_2 = 0;
  for(int i = 0; i < n; i++){
    gammahat_j[i] = data[i];
    Gammahat_j[i] = data[n + i];
    sigma2_X[i] = pow(data[2*n + i], 2);
    sigma2_Y[i] = pow(data[3*n + i], 2);
    psi2_j[i] = sigma2_X[i] + hyper_gammaj_psi2;
    tau2_j[i] = sigma2_Y[i] + rhyper_Gammaj_tau2;
    h2_j[i] = (pow(beta_p, 2))*hyper_gammaj_psi2 + tau2_j[i];
    sum_A_beta_1 += 1/psi2_j[i];
    sum_A_beta_2 += 1/h2_j[i];
    sum_B_beta_1 += gammahat_j[i]/psi2_j[i];
    sum_B_beta_2 += Gammahat_j[i]/(beta_p*h2_j[i]);
  }

  GetRNGstate();

  double A_beta = 0, B_beta = 0;
  double prob = 0, ran_unif = 0;
  int accpeted_beta = 0;
  double* lpost_beta_prop = new double;
  double* lpost_beta_old = new double;
  accept[0] = 0;
  beta_old = beta_p;
  while (niter < totiter){
    niter++;
    // generate gamma by Gibbs sampling
    A_beta = sum_A_beta_1 + pow(beta_old, 2)*sum_A_beta_2 + 1/rhyper_gamma_var;
    B_beta = sum_B_beta_1 + pow(beta_old, 2)*sum_B_beta_2 +
      rhyper_gamma_mean/rhyper_gamma_var;
    gamma_chain[(niter - 1)] = R::rnorm(B_beta/A_beta, sqrt(1/A_beta));
    
    // generate beta by random walk Metropolis-Hastings
    // beta_prop = R::runif(beta_old - sigma2_beta, beta_old + sigma2_beta);  // uniform proposal
    beta_prop = R::rnorm(beta_old, sqrt(sigma2_beta));  // normal proposal
    logpost_beta(lpost_beta_prop, beta_prop, gamma_chain[(niter - 1)], rhyper_beta_mean,
      rhyper_beta_var, hyper_gammaj_psi2, rhyper_Gammaj_tau2, n, Gammahat_j, sigma2_Y);
    logpost_beta(lpost_beta_old, beta_old, gamma_chain[(niter - 1)], rhyper_beta_mean,
      rhyper_beta_var, hyper_gammaj_psi2, rhyper_Gammaj_tau2, n, Gammahat_j, sigma2_Y);
    prob = exp(*lpost_beta_prop - *lpost_beta_old);
    
    ran_unif = R::runif(0, 1);
    accpeted_beta = (ran_unif < prob) ? 1 : 0;
    accept[0] += accpeted_beta;
    beta_old = (accpeted_beta == 1) ? beta_prop : beta_old;
    beta_chain[(niter - 1)] = beta_old;
    accept_beta = accept[0];  // [SV: this is redundant --> simplify in the future]

    // print the information
    if(((niter % 500) == 0) && verbose){
      REprintf("   iteration %d/%d ==> acceptance beta: %1.4f\n", niter, totiter, accept_beta/(G*niter));
    }

    // calculate the loglikelihood, logprior and logposterior for the sampled parameter values
    logprior[niter - 1] = 0;
    // for(int g = 0; g < G; g++){
    //   for(int j = 0; j < p; j++){
    //     for(int i = 0; i < n; i++){
    //       z_g[i + n*j] = z[i + n*j + n*p*g];
    //     }
    //   }
    //   for(int j = 0; j < p; j++){
    //     sigma_g[j + p*j] = eta[g];
    //   }
    //   dmultinorm(lprior_z, z_g, mean_g, sigma_g, n, p, 1);
    //   sum_lprior_z = 0;
    //   for(int i = 0; i < n; i++){
    //     sum_lprior_z += lprior_z[i];
    //   }
    //   lprior_alpha = R::dnorm(alpha[g], 0, sqrt(sigma2[g]), 1);
    //   sum_z2 = 0;
    //   for(int j = 0; j < p; j++){
    //     for(int i = 0; i < n; i++){
    //       sum_z2 += z_g[i + n*j]*z_g[i + n*j];
    //     }
    //   }
    //   dinvgamma(lprior_eta, &eta[g], hyper_eta_a[g], hyper_eta_b[g], 1, 1);
    //   dinvgamma(lprior_sigma2, &sigma2[g], hyper_sigma2_a, hyper_sigma2_b, 1, 1);
    //   logprior[niter - 1] += sum_lprior_z + lprior_alpha + *lprior_eta + *lprior_sigma2;
    // }
    // ddirichlet(lprior_lambda, lambda, hyper_lambda, 1, G, 1);
    // logprior[niter - 1] += *lprior_lambda;
    // loglik_bayesmr(&loglik[niter - 1], Dm, z, alpha, d_sigma2, lambda, x, n, p, S, G, "binomial");
    loglik[niter - 1] = 0;
    logpost[niter - 1] = loglik[niter - 1] + logprior[niter - 1];
    R_CheckUserInterrupt();
  }

  PutRNGstate();

  delete[] gammahat_j;
  delete[] Gammahat_j;
  delete[] sigma2_X;
  delete[] sigma2_Y;
  delete[] psi2_j;
  delete[] tau2_j;
  delete[] h2_j;

  delete   lpost_beta_prop;
  delete   lpost_beta_old;
}
