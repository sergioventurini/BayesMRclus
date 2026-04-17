library(BayesMRclus)

# prepare data
data("hdl_chd", package = "BayesMRclus")

dat <- subset(hdl_chd, pval_selection < 5e-8)

n <- nrow(dat)
hdl_chd <- new("bayesmr_data", data = dat, n = n, reorientation = TRUE)
dat_tmp <- hdl_chd@data
# summary(hdl_chd)
# plot(hdl_chd, se = TRUE)

# simulation setup
prm.prop <- list(beta = sqrt(0.1), psi = sqrt(.15), tau = sqrt(1))
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
alpha_tau <- 0  # no pleiotropy
prior <- bayesmr_prior(psi = list(alpha = alpha_psi, nu = 3),
                       tau = list(alpha = alpha_tau, nu = 3),
                       gamma = list(mean = 0, var = 100),
                       beta = list(mean = 0, var = 100),
                       alpha = list(a = 2, b = 1))

# MCMC simulation
res_BNPM_MR <- bayesmr_mix_het(hdl_chd, control = control, prior = prior)

# Results ======================================================================

dat_post <- hdl_chd@data
dat_post[, 3] <- dat_post[, 3]^2
dat_post[, 4] <- dat_post[, 4]^2
names(dat_post) <- c("gamma_hat", "Gamma_hat", "sigma2X", "sigma2Y")

results <- run_bnpmr_postprocess(
  res       = res_BNPM_MR,
  burnin    = burnin/control$thin,
  dat       = dat_post,  # needed for pleiotropy_check; NULL to skip
  fig_dir   = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/figures/MR",  # NULL to suppress figure output
  ci        = 0.95,
  ordering  = "post_mean",
  loss      = "binder",
  verbose   = FALSE,
  snp_names = rownames(dat_post),
  name_cex  = 0.5,
  pleiotropy = (prior$tau$alpha > 0)
)

# results$K_posterior          # posterior distribution of K
# results$binder$partition     # representative partition (1-indexed)
# results$clusters             # cluster-specific effect summaries
# results$bma                  # per-SNP BMA posterior effects
# results$pleiotropy$mean_T    # dispersion statistic

# mcmc_sim <- cbind(as.data.frame(results$out[1:3]), K = apply(results$out$xi, 1, max))
# write.csv(mcmc_sim, file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/mcmc_sim.csv")
# write.csv(results$K_posterior , file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/Kpost.csv")
# write.csv(results$binder$partition, file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/binder.csv")
# write.csv(results$clusters, file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/clus.csv")
# write.csv(results$bma, file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/bma.csv")

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
cl_labels <- c("1"="Negative", "2"="Positive", "3"="Null")

cl_agr_res <- compare_partitions(
  part_mrclust           = cluster_agr$MRCLUST,
  part_bnpmr             = cluster_agr$BNPM_MR,
  snp_names              = mrclust_res$cluster$SNP,
  wald                   = wald,
  se_wald                = se_w,
  cluster_labels_mrclust = cl_labels,
  fig_dir                = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/figures/MR",
  csv_dir                = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR",
  verbose                = FALSE
)

# cl_agr_res$agreement$measures            # ARI, NVI, FMI
# cl_agr_res$agreement$contingency$table   # raw contingency table
# cl_agr_res$snp_table                     # per-SNP concordance

save.image(file = "/Users/Sergio/dev/BayesMRclus/test/examples_mix_het/example_hdl_chd/output/MR/out_MR.rda")

# MR-PATH ======================================================================

library(MRPATH)  # it seems to work only with at most K = 2 clusters

data_mrpath <- hdl_chd@data
data_mrpath <- cbind(SNP = rownames(data_mrpath), data_mrpath)
colnames(data_mrpath) <- gsub("_", ".", colnames(data_mrpath))

# model selection
K_max <- 2

fit <- vector(mode = "list", length = K_max)
BIC <- numeric(K_max)
for (K in 1:K_max) {
  InitVals <- MRPATH_optimizeInitVals(K, data_mrpath, Nreps = 10)$initVals
  fit[[K]] <- MR_PATH(K, data_mrpath, InitVals)
  BIC[K] <- MRPATH_BIC(data_mrpath, K, fit[[K]])
}
K_best <- which.min(BIC)
fit_best <- fit[[K_best]]

### Compute SNP-specific cluster membership probabilities 
clustermemb_prob <- computeClusterMembProb(data_mrpath, MCEM_fit = fit_best)
clustermemb <- apply(clustermemb_prob, 1, which.max)
MRPATH_barplot(data_mrpath, fit_best, ret.snps = FALSE)

# as in the paper by Iong et al. (2024) ==> in the MRPATH examples
# they DO NOT reorient the SNPs!
colnames(dat) <- gsub("_", ".", colnames(dat))
K <- 2
initVals <- list("m_X" = mean(dat$beta.exposure),
                 "lambdaX" = sd(dat$beta.exposure),
                 "pis" = rep(1/K, K),
                 "mus" = c(-.9, .4),
                 "sds" = c(.9, .4))
set.seed(1406)
MCEM_fit <- MR_PATH(K, dat, initVals)
clustermemb_prob <- computeClusterMembProb(dat, MCEM_fit = MCEM_fit)
clustermemb <- apply(clustermemb_prob, 1, which.max)
MRPATH_barplot(data_mrpath, MCEM_fit, ret.snps = FALSE)

# BNPM-MR vs. MR-Path
table(results$binder$partition, clustermemb)
