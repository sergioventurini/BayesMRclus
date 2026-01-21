library(BayesMRclus)

# prepare data
data("hdl_chd", package = "BayesMRclus")
hdl_chd <- subset(hdl_chd, pval.selection < 5e-8)
data <- data.frame(beta_exposure = hdl_chd[, "beta.exposure"],
                   beta_outcome = hdl_chd[, "beta.outcome"],
                   se_exposure = hdl_chd[, "se.exposure"],
                   se_outcome = hdl_chd[, "se.outcome"])

n <- nrow(data)
iongdata <- new("bayesmr_data", data = data, n = n, harmonization = TRUE)
# summary(iongdata)
# plot(iongdata, se = TRUE)

# simulation setup
prm.prop <- list(beta = .2)
burnin <- 10000
nsim <- 20000

# seed <- 2301
# set.seed(seed)

nchains <- 3
control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                beta.m = 2, nchains = nchains, thin = 1,
                random.start = TRUE, verbose = TRUE,
                store.burnin = TRUE, threads = ifelse(
                  nchains <= parallel::detectCores(),
                  nchains, parallel::detectCores() - 1),
                parallel = "snow")

# tau2_est <- tau2_dl(data_tmp, secondorder = TRUE)
prior <- bayesmr_prior(gammaj = list(psi2 = .0001),
                       Gammaj = list(tau2 = .0001),
                       gamma = list(mean = 0, var = 1e1),
                       beta = list(mean = 0, var = 1e1),
                       alpha = list(a = 2, b = 2))

# MCMC simulation
res_BMR <- bayesmr_mix(iongdata, control = control, prior = prior)

plot(res_BMR@results[[1]]@gamma.chain, type = "l")
lines(res_BMR@results[[2]]@gamma.chain, type = "l", col = gray(.6))
lines(res_BMR@results[[3]]@gamma.chain, type = "l", col = gray(.9))
plot(res_BMR@results[[1]]@alpha.chain, type = "l")
lines(res_BMR@results[[2]]@alpha.chain, type = "l", col = gray(.6))
lines(res_BMR@results[[3]]@alpha.chain, type = "l", col = gray(.9))
str(res_BMR@results[[1]]@beta.chain)
str(res_BMR@results[[1]]@xi.chain)

# MCMC summaries/graphs
summary(res_BMR)
plot(res_BMR, what = "trace", regex_pars = "beta")
plot(res_BMR, what = "trace", regex_pars = "gamma")
# plot(res_BMR, what = "hist_by_chain", regex_pars = "beta")
# plot(res_BMR, what = "dens_chains", regex_pars = "beta")
# plot(res_BMR, what = "intervals", regex_pars = "beta")
# plot(res_BMR, what = "acf_bar", regex_pars = "beta")

# ch <- 1
# beta_c <- res_BMR@results[[ch]]@beta.chain[1001:3000]
# plot(beta_c, type = "l")
# gamma_c <- res_BMR@results[[ch]]@gamma.chain[1001:3000]
# plot(gamma_c, type = "l")
# beta_gamma_c <- beta_c*gamma_c
# plot(beta_gamma_c, type = "l")
# summary(beta_gamma_c)
# plot(beta_c, gamma_c, pch = 20)

# coda analyses
library(coda)
res_BMR_sub <- subset(res_BMR, regex_pars = c("gamma", "beta"))
plot(res_BMR_sub)
cumuplot(res_BMR_sub)
gelman.plot(res_BMR_sub)
geweke.plot(res_BMR_sub)
raftery.diag(res_BMR_sub)
HPDinterval(res_BMR_sub)
heidel.diag(res_BMR_sub)
densplot(res_BMR_sub)

# post-processing
post_sim <- function(xi_post, burnin) {
  p <- ncol(xi_post)
  simil_mat <- matrix(0, nrow = p, ncol = p)
  
  S <- nrow(xi_post)

  for(s in (burnin + 1):S){
    simil_mat <- simil_mat + (matrix(xi_post[s, ], nrow = p, ncol = p) == t(matrix(xi_post[s, ], nrow = p, ncol = p)))*1
  }
  
  simil_probs <-simil_mat/(S - burnin)
  
  return(simil_probs)
}

psm1 <- post_sim(res_BMR@results[[1]]@xi.chain, burnin = burnin)
psm2 <- post_sim(res_BMR@results[[2]]@xi.chain, burnin = burnin)
psm3 <- post_sim(res_BMR@results[[3]]@xi.chain, burnin = burnin)
plot(psm1, psm2) ## good agreement
plot(psm1, psm3) ## good agreement
plot(psm2, psm3) ## good agreement

library(mcclust.ext)
library(mclust)
library(fields)

## Wade Gaharamani BA 2018

vi =  minVI(psm1)$cl
table(vi)

ordered_ind_vi <- unlist(sapply(length(table(vi)):1,
                                function(k) c(which(vi == k))))

simil_ordered_vi <- psm1[ordered_ind_vi, ordered_ind_vi]

colori <- colorRampPalette(c('white','black'))

p <- ncol(simil_ordered_vi)

par(mar = c(5,5,2,2), mgp = c(3, 1, 0))
image.plot(simil_ordered_vi, col = colori(100),
           zlim = c(0,1), cex.sub = 3, axes = FALSE,
           horizontal = FALSE, legend.shrink = 1,
           xlab = "SNP (j)", ylab = "SNP (j')",
           cex.axis = 4 )
axis(1, at = seq(0, 1, length = p), lab = 1:p, las = 2)
axis(2, at = seq(0, 1, length = p), lab = 1:p, las = 2)

plot(sapply(1:nrow(res_BMR@results[[1]]@xi.chain),
            function(s) length(unique(res_BMR@results[[1]]@xi.chain[s, ]))),
     ylab = "Number of clusters", type = "l", xlab = "MCMC iteration")

S <- nrow(res_BMR@results[[1]]@beta.chain)

# Recover SNP-specific causal effects

beta_snps <- sapply(1:S, function(s) res_BMR@results[[1]]@beta.chain[s, ][res_BMR@results[[1]]@xi.chain[s, ]])
boxplot(t(beta_snps[, burnin:S]), outline = FALSE)

# arranged according to the estimated partition

boxplot(t(beta_snps[ordered_ind_vi, burnin:S]), outline = FALSE,
        col = c(rep("lightblue", sum(vi == 1)), rep("pink", sum(vi == 2))))

# Sorted in a more meaningful way

means <- colMeans(t(beta_snps[, burnin:S]))
ord <- order(means)

hdl_chd$SNP[ord]
plot(means[ord])

par(mar = c(8, 8, 2, 2), mgp = c(6, 1, 0))
image.plot(psm1[ord, ord], col = colori(100), zlim = c(0,1), cex.sub = 3, 
           axes = FALSE, horizontal = FALSE, legend.shrink = 1,
           xlab = "SNP (j)", ylab = "SNP (j')", cex.axis = 4 )
axis(1, at = seq(0, 1, length = p), lab = hdl_chd$SNP[ord], las = 2)
axis(2, at = seq(0, 1, length = p), lab = hdl_chd$SNP[ord], las = 2)

###

cols <- ifelse(vi == 1, "lightblue", "pink")
beta_post_ord <- t(beta_snps[ord, burnin:S])
colnames(beta_post_ord) <- hdl_chd$SNP[ord]
par(mar = c(8, 5, 2, 2))
boxplot(beta_post_ord, outline = FALSE, col = "gray", las = 2,
        ylim = c(-1.5, 1.5),
        xlab = "",
        ylab = expression(beta[j]))
title(xlab = "SNP (j)", line = 6)

# abline(h = beta_no_mixture, col = "black", lty = 2)
abline(h = 0, col = "black", lty = 2)

###

par(mar = c(8, 5, 2, 2))
boxplot(beta_post_ord, outline = FALSE, col = cols[ord], las = 2,
        ylim = c(-1.5, 1.5),
        xlab = "",
        ylab = expression(beta[j]))
title(xlab = "SNP (j)", line = 6)

# abline(h = beta_no_mixture, col = "black", lty = 2)
abline(h = 0, col = "black", lty = 2)
