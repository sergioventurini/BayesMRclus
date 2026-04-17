# BNPM-MR =====================================================================

library(BayesMRclus)

# prepare data
data("whradjbmi_insulin", package = "BayesMRclus")

dat <- data.frame(SNP = whradjbmi_insulin[, "SNP"],
                  beta_exposure = whradjbmi_insulin[, "beta.exposure"],
                  beta_outcome = whradjbmi_insulin[, "beta.outcome"],
                  se_exposure = whradjbmi_insulin[, "se.exposure"],
                  se_outcome = whradjbmi_insulin[, "se.outcome"])

n <- nrow(whradjbmi_insulin)
whradjbmi_insulin <- new("bayesmr_data", data = dat, n = n, reorientation = TRUE)
dat_tmp <- whradjbmi_insulin@data
# summary(whradjbmi_insulin)
# plot(whradjbmi_insulin, se = TRUE)

# simulation setup
prm.prop <- list(beta = sqrt(0.1), psi = sqrt(.025), tau = sqrt(1))
burnin <- 20000
nsim <- 100000

nchains <- 3
seed <- 1103 #round(runif(1)*10000)
control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                psi.prop = prm.prop[["psi"]], tau.prop = prm.prop[["tau"]],
                beta.m = 2, nchains = nchains, thin = 40,
                seed = seed, random_start = TRUE, verbose = TRUE,
                store.burnin = TRUE, threads = ifelse(
                  nchains <= parallel::detectCores() - 1,
                  nchains, parallel::detectCores() - 1),
                parallel = "snow")

alpha_psi <- median(dat_tmp[, "se_exposure"])
alpha_tau <- 0
prior <- bayesmr_prior(psi = list(alpha = alpha_psi, nu = 3),
                       tau = list(alpha = alpha_tau, nu = 3),
                       gamma = list(mean = 0, var = 10),
                       beta = list(mean = 0, var = 150),
                       alpha = list(a = 3, b = 1))

# MCMC simulation
res_BNPM_MR <- bayesmr_mix_het(whradjbmi_insulin, control = control, prior = prior)

# Results ======================================================================

dat_post <- whradjbmi_insulin@data
dat_post[, 3] <- dat_post[, 3]^2
dat_post[, 4] <- dat_post[, 4]^2
names(dat_post) <- c("gamma_hat", "Gamma_hat", "sigma2X", "sigma2Y")

results <- run_bnpmr_postprocess(
  res       = res_BNPM_MR,
  burnin    = burnin/control$thin,
  dat       = dat_post,  # needed for pleiotropy_check; NULL to skip
  fig_dir   = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_whradjbmi_insulin/figures/MR",  # NULL to suppress figure output
  ci        = 0.95,
  ordering  = "post_mean",
  loss      = "binder",
  verbose   = FALSE,
  snp_names = NULL, #rownames(dat_post),
  name_cex  = 0.3,
  pleiotropy = (prior$tau$alpha > 0)
)

# results$K_posterior          # posterior distribution of K
# results$binder$partition     # representative partition (1-indexed)
# results$clusters             # cluster-specific effect summaries
# results$bma                  # per-SNP BMA posterior effects
# results$pleiotropy$mean_T    # dispersion statistic

# MR-CLUST =====================================================================

if (!exists("mrclust_res")) {
  library(ggplot2); library(dplyr); library(mrclust)
  
  seed <- 3110
  set.seed(seed)
  
  dat_mrclust <- dat_tmp
  colnames(dat_mrclust) <- gsub("_", ".", colnames(dat_mrclust))
  
  mrclust_res <- mrclust_results(dat_mrclust)
  # mrclust_res$p1
  mrclust_res$tab
  mrclust_res$cluster_class
}

# Cluster agreement ============================================================

cluster_agr <- data.frame(SNP = dat$SNP,
                          BNPM_MR = results$binder$partition,
                          MRCLUST = mrclust_res$cluster$cluster,
                          MRCLUST_class = mrclust_res$cluster$cluster_class)
table(cluster_agr$BNPM_MR, cluster_agr$MRCLUST)

# Wald ratios and their SEs
wald  <- dat_post$Gamma_hat / dat_post$gamma_hat
se_w  <- abs(sqrt(dat_post$sigma2Y) / dat_post$gamma_hat)

#  # Descriptive labels matching the MRClust integer codes above
cl_labels <- c("1"="Negative", "2"="Positive", "3"="Null", "4"="Junk")

cl_agr_res <- compare_partitions(
  part_mrclust           = cluster_agr$MRCLUST,
  part_bnpmr             = cluster_agr$BNPM_MR,
  snp_names              = NULL, #mrclust_res$cluster$SNP,
  wald                   = wald,
  se_wald                = se_w,
  cluster_labels_mrclust = cl_labels,
  fig_dir                = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_whradjbmi_insulin/figures/MR",
  csv_dir                = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_whradjbmi_insulin/output/MR",
  verbose                = FALSE
)

# cl_agr_res$agreement$measures            # ARI, NVI, FMI
# cl_agr_res$agreement$contingency$table   # raw contingency table
# cl_agr_res$snp_table                     # per-SNP concordance

# MR-PATH ======================================================================

library(MRPATH)  # it seems to work only with at most K = 2 clusters

data_mrpath <- whradjbmi_insulin@data
colnames(data_mrpath) <- gsub("_", ".", colnames(data_mrpath))

# model selection
K_max <- 2
modSel <- MRPATH_selectModel(data_mrpath, K_range = 1:K_max, Nreps = 20,
                             altModel = FALSE, verbose = TRUE)
K <- modSel$bestK

# initial values
initVals <- MRPATH_optimizeInitVals(K, data_mrpath, Nreps = 10, verbose = FALSE,
                                    altModel = FALSE, init_seed = 3110)$initVals
# initVals <- list("m_X" = mean(data_mrpath$beta.exposure),
#                  "lambdaX" = sd(data_mrpath$beta.exposure),
#                  "pis" = rep(1/K, K),
#                  "mus" = runif(K, -2, 2),
#                  "sds" = runif(K, 0, 1))

### Run MC-EM algorithm
fit <- MR_PATH(K, data_mrpath, initVals)

### Compute SNP-specific cluster membership probabilities 
clustermemb_prob <- computeClusterMembProb(data_mrpath, MCEM_fit = fit)
# MRPATH_barplot(data_mrpath, fit, ret.snps = FALSE)
clustermemb <- apply(clustermemb_prob, 1, which.max)
