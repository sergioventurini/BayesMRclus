library(LaplacesDemon)
library(MCMCpack)
library(BayesMRclus)

# set.seed(101)
nsim <- 500
p <- 100
K <- 2
pi_k <- c(.5, .5)
sum(pi_k) == 1
mu_k <- c(-.5, .5)
beta_i <- sample(mu_k, size = p, replace = TRUE, prob = pi_k)
table(beta_i)
lambda_x <- 10/sqrt(p)
theta_x <- rnorm(p, 0, lambda_x)
sigma2_x <- MCMCpack::rinvgamma(p, shape = 9, scale = 0.002)
sigma2_y <- MCMCpack::rinvgamma(p, shape = 9, scale = 0.002)
tau_0 <- (2/p)*sum(sqrt(sigma2_y))

# scenario 1
alpha_i <- rnorm(p, 0, tau_0)

# scenario 2
# alpha_i <- NA

# scenario 3
# alpha_i <- NA

theta_hat_x <- rnorm(p, theta_x, sqrt(sigma2_x))
theta_hat_y <- rnorm(p, alpha_i + beta_i*theta_x, sqrt(sigma2_y))

data <- data.frame(beta_exposure = theta_hat_x,
                   beta_outcome = theta_hat_y,
                   se_exposure = sqrt(sigma2_x),
                   se_outcome = sqrt(sigma2_y))
simdata <- new("bayesmr_data", data = data, n = nrow(data), reorientation = TRUE)
# summary(simdata)
# plot(simdata, se = TRUE)
hist(data$beta_outcome/data$beta_exposure, breaks = 20)

###

rlaplace <- function(n, mu = 0, b = 1) {
  u <- runif(n, -0.5, 0.5)
  return(mu - b * sign(u) * log(1 - 2*abs(u)))
}

rinvgamma <- function(n, shape, scale = 1) {
  return(1 / rgamma(n, shape = shape, rate = scale))
}

simulate_MRPath <- function(
    p       = 100,
    K       = 2,
    pi      = c(0.5, 0.5),
    mu      = c(-0.5, 0.5),
    lambda  = 10/sqrt(p),
    alpha_type = c("normal","laplace","idiosyncratic")
) {
  alpha_type <- match.arg(alpha_type)
  
  ## cluster membership
  z <- sample(1:K, size = p, replace = TRUE, prob = pi)
  beta <- mu[z]
  
  ## θ_Xi ∼ N(0, λ_x^2)
  theta_X <- rnorm(p, mean = 0, sd = lambda)
  
  ## measurement errors σ_X^2, σ_Y^2 ∼ InvGamma(9, .002)
  sigma2_X <- rinvgamma(p, shape = 9, scale = .002)
  sigma2_Y <- rinvgamma(p, shape = 9, scale = .002)
  
  ## τ0 = (2/p) * sum σ_Y
  tau0 <- (2/p) * sum(sqrt(sigma2_Y))
  
  ## pleiotropy α_i
  if (alpha_type == "normal") {
    alpha <- rnorm(p, mean = 0, sd = tau0)
  } else if (alpha_type == "laplace") {
    alpha <- tau0 * rlaplace(p, b = 1)
  } else if (alpha_type == "idiosyncratic") {
    alpha <- rnorm(p, mean = 0, sd = tau0)
    idx <- sample(1:p, size = round(0.10*p))   # 10% spikes
    alpha[idx] <- rnorm(length(idx), mean = 5*tau0, sd = tau0)
  }
  
  ## generate observed summary stats:
  ##  (θ̂_X , θ̂_Y)
  theta_hat_X <- rnorm(p, mean = theta_X, sd = sqrt(sigma2_X))
  theta_hat_Y <- rnorm(p, mean = alpha + beta*theta_X, sd = sqrt(sigma2_Y))
  
  return(list(
    data = data.frame(beta_exposure = theta_hat_X,
                      beta_outcome = theta_hat_Y,
                      se_exposure = sqrt(sigma2_X),
                      se_outcome = sqrt(sigma2_Y)),
    beta = beta,
    z = z,
    theta_X = theta_X,
    alpha = alpha,
    se_exposure = sqrt(sigma2_X),
    se_outcome = sqrt(sigma2_Y),
    tau0 = tau0
  ))
}

# set.seed(101)
data <- simulate_MRPath(alpha_type = "normal")$data
simdata <- new("bayesmr_data", data = data, n = nrow(data), reorientation = TRUE)
# summary(simdata)
# plot(simdata, se = TRUE)
hist(data$beta_outcome/data$beta_exposure, breaks = 20)
