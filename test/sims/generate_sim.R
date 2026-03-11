# correctly specified: moderate heterogeneity, mild pleiotropy
simulate_scenario1 <- function(p = 150, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # True hyperparameters
  gamma_true <- 0.08
  psi_true   <- 0.03
  tau_true   <- 0.02
  
  # Measurement SEs (fixed for simplicity)
  s_gamma_hat <- rep(0.01, p)
  s_Gamma_hat <- rep(0.02, p)
  
  # Cluster allocation
  z <- rbinom(p, 1, 0.6)
  beta_true <- ifelse(z == 1, 0.25, 0.00)
  
  # Generate latent exposure effects
  gamma_j <- rnorm(p, gamma_true, psi_true)
  
  # Generate summary statistics
  gamma_hat <- rnorm(p, gamma_j, s_gamma_hat)
  Gamma_hat <- rnorm(p, beta_true * gamma_j, 
                     sqrt(s_Gamma_hat^2 + tau_true^2))
  
  data.frame(
    gamma_hat = gamma_hat,
    Gamma_hat = Gamma_hat,
    s_gamma_hat = s_gamma_hat,
    s_Gamma_hat = s_Gamma_hat,
    beta_true = beta_true
  )
}

# correctly specified: high complexity --> three mechanisms + strong pleiotropy
simulate_scenario2 <- function(p = 200, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  gamma_true <- 0.07
  psi_true   <- 0.05
  tau_true   <- 0.05
  
  # Heterogeneous SEs
  s_gamma_hat <- runif(p, 0.008, 0.015)
  s_Gamma_hat <- runif(p, 0.015, 0.03)
  
  # Cluster allocation
  z <- sample(1:3, p, replace = TRUE, 
              prob = c(0.5, 0.3, 0.2))
  
  beta_values <- c(0.30, 0.10, -0.15)
  beta_true <- beta_values[z]
  
  # Latent exposure effects
  gamma_j <- rnorm(p, gamma_true, psi_true)
  
  # Summary statistics
  gamma_hat <- rnorm(p, gamma_j, s_gamma_hat)
  Gamma_hat <- rnorm(p, beta_true * gamma_j,
                     sqrt(s_Gamma_hat^2 + tau_true^2))
  
  data.frame(
    gamma_hat = gamma_hat,
    Gamma_hat = Gamma_hat,
    s_gamma_hat = s_gamma_hat,
    s_Gamma_hat = s_Gamma_hat,
    beta_true = beta_true
  )
}

# misspecified: latent mechanistic pleiotropy
simulate_scenario3 <- function(p = 200, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # True hyperparameters
  gamma_true <- 0.08
  psi_true   <- 0.04
  sigma_eps  <- 0.02
  
  # Measurement SEs
  s_gamma_hat <- runif(p, 0.008, 0.015)
  s_Gamma_hat <- runif(p, 0.015, 0.03)
  
  # Cluster allocation
  z <- sample(1:3, p, replace = TRUE,
              prob = c(0.5, 0.3, 0.2))
  
  beta_vals    <- c(0.30, 0.05, -0.20)
  lambda_vals  <- c(0.10, -0.08, 0.12)
  
  beta_true   <- beta_vals[z]
  lambda_true <- lambda_vals[z]
  
  # Latent pathway activity
  U_j <- rnorm(p, 0, 1)
  
  # Exposure effects
  gamma_j <- rnorm(p, gamma_true, psi_true)
  
  # True outcome association (misspecified relative to BNPM)
  Gamma_j_true <- beta_true * gamma_j + 
                  lambda_true * U_j + 
                  rnorm(p, 0, sigma_eps)
  
  # Summary statistics
  gamma_hat <- rnorm(p, gamma_j, s_gamma_hat)
  Gamma_hat <- rnorm(p, Gamma_j_true, s_Gamma_hat)
  
  data.frame(
    gamma_hat = gamma_hat,
    Gamma_hat = Gamma_hat,
    s_gamma_hat = s_gamma_hat,
    s_Gamma_hat = s_Gamma_hat,
    beta_true = beta_true,
    cluster_true = z
  )
}

# misspecified: latent pathways + weak instruments
simulate_scenario4 <- function(p = 250, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Proportion weak
  weak_prob <- 0.4
  
  # Cluster allocation (mechanistic pathways)
  z <- sample(1:3, p, replace = TRUE,
              prob = c(0.5, 0.3, 0.2))
  
  beta_vals   <- c(0.30, 0.05, -0.20)
  lambda_vals <- c(0.10, -0.08, 0.12)
  
  beta_true   <- beta_vals[z]
  lambda_true <- lambda_vals[z]
  
  # Weak vs strong instruments
  weak_flag <- rbinom(p, 1, weak_prob)
  
  gamma_j <- ifelse(
    weak_flag == 1,
    rnorm(p, 0.01, 0.02),   # weak
    rnorm(p, 0.08, 0.04)    # strong
  )
  
  # Latent pathway factor
  U_j <- rnorm(p, 0, 1)
  
  # True structural outcome
  Gamma_j_true <- beta_true * gamma_j +
                  lambda_true * U_j +
                  rnorm(p, 0, 0.02)
  
  # Heterogeneous SEs
  s_gamma_hat <- runif(p, 0.01, 0.02)
  s_Gamma_hat <- runif(p, 0.015, 0.03)
  
  gamma_hat <- rnorm(p, gamma_j, s_gamma_hat)
  Gamma_hat <- rnorm(p, Gamma_j_true, s_Gamma_hat)
  
  data.frame(
    gamma_hat = gamma_hat,
    Gamma_hat = Gamma_hat,
    s_gamma_hat = s_gamma_hat,
    s_Gamma_hat = s_Gamma_hat,
    beta_true = beta_true,
    cluster_true = z,
    weak_flag = weak_flag
  )
}

###

data_sc1 <- simulate_scenario1(p = 150, seed = 101)
dat <- data.frame(beta_exposure = data_sc1[, "gamma_hat"],
                  beta_outcome = data_sc1[, "Gamma_hat"],
                  se_exposure = data_sc1[, "s_gamma_hat"],
                  se_outcome = data_sc1[, "s_Gamma_hat"])

dat_sc1 <- new("bayesmr_data", data = dat, n = nrow(dat), harmonization = TRUE)
# summary(dat_sc1)
# plot(dat_sc1, se = TRUE)
hist(dat_sc1@data$beta_outcome/dat_sc1@data$beta_exposure, breaks = 20)

# library(MendelianRandomization)
# mr_ivw(mr_input(dat_sc1@data$beta_exposure, dat_sc1@data$se_exposure,
#                 dat_sc1@data$beta_outcome, dat_sc1@data$se_outcome))
# mr_allmethods(mr_input(dat_sc1@data$beta_exposure, dat_sc1@data$se_exposure,
#                        dat_sc1@data$beta_outcome, dat_sc1@data$se_outcome))

#

data_sc2 <- simulate_scenario2(p = 150, seed = 201)
dat <- data.frame(beta_exposure = data_sc2[, "gamma_hat"],
                  beta_outcome = data_sc2[, "Gamma_hat"],
                  se_exposure = data_sc2[, "s_gamma_hat"],
                  se_outcome = data_sc2[, "s_Gamma_hat"])

dat_sc2 <- new("bayesmr_data", data = dat, n = nrow(dat), harmonization = TRUE)
# summary(dat_sc2)
# plot(dat_sc2, se = TRUE)
hist(dat_sc2@data$beta_outcome/dat_sc2@data$beta_exposure, breaks = 20)

# library(MendelianRandomization)
# mr_ivw(mr_input(dat_sc2@data$beta_exposure, dat_sc2@data$se_exposure,
#                 dat_sc2@data$beta_outcome, dat_sc2@data$se_outcome))
# mr_allmethods(mr_input(dat_sc2@data$beta_exposure, dat_sc2@data$se_exposure,
#                        dat_sc2@data$beta_outcome, dat_sc2@data$se_outcome))

#

data_sc3 <- simulate_scenario3(p = 150, seed = 301)
dat <- data.frame(beta_exposure = data_sc3[, "gamma_hat"],
                  beta_outcome = data_sc3[, "Gamma_hat"],
                  se_exposure = data_sc3[, "s_gamma_hat"],
                  se_outcome = data_sc3[, "s_Gamma_hat"])

dat_sc3 <- new("bayesmr_data", data = dat, n = nrow(dat), harmonization = TRUE)
# summary(dat_sc3)
# plot(dat_sc3, se = TRUE)
abline(a = 0, b = 0.3)
hist(dat_sc3@data$beta_outcome/dat_sc3@data$beta_exposure, breaks = 20)

#

data_sc4 <- simulate_scenario4(p = 150, seed = 401)
dat <- data.frame(beta_exposure = data_sc4[, "gamma_hat"],
                  beta_outcome = data_sc4[, "Gamma_hat"],
                  se_exposure = data_sc4[, "s_gamma_hat"],
                  se_outcome = data_sc4[, "s_Gamma_hat"])

dat_sc4 <- new("bayesmr_data", data = dat, n = nrow(dat), harmonization = F)
# summary(dat_sc4)
# plot(dat_sc4, se = TRUE)
hist(dat_sc4@data$beta_outcome/dat_sc4@data$beta_exposure, breaks = 20)
