library(BayesMRclus)

# prepare data
data("hdl_cad", package = "BayesMRclus")

dat <- data.frame(SNP = hdl_cad[, "SNP"],
                  beta_exposure = hdl_cad[, "beta.exposure"],
                  beta_outcome = hdl_cad[, "beta.outcome"],
                  se_exposure = hdl_cad[, "se.exposure"],
                  se_outcome = hdl_cad[, "se.outcome"])

n <- nrow(hdl_cad)
hdlcad_data <- new("bayesmr_data", data = dat, n = n, harmonization = TRUE)
data_tmp <- hdlcad_data@data
# summary(hdlcad_data)
# plot(hdlcad_data, se = TRUE)

# simulation setup
prm.prop <- list(beta = 1.5, psi = .035, tau = .035)
burnin <- 10000
nsim <- 20000

# seed <- 2301
# set.seed(seed)

nchains <- 3
control <- list(burnin = burnin, nsim = nsim, beta.prop = prm.prop[["beta"]],
                psi.prop = prm.prop[["psi"]], tau.prop = prm.prop[["tau"]],
                beta.m = 2, nchains = nchains, thin = 5,
                random_start = TRUE, verbose = TRUE,
                store.burnin = FALSE, threads = ifelse(
                  nchains <= parallel::detectCores() - 1,
                  nchains, parallel::detectCores() - 1),
                parallel = "snow")

alpha_psi <- 0.01 #0.05 #alpha_eb(se = data_tmp[, "se_exposure"])
alpha_tau <- 2*median(data_tmp[, "se_outcome"]) #2*alpha_psi #alpha_eb(se = data_tmp[, "se_outcome"])
prior <- bayesmr_prior(psi = list(alpha = alpha_psi, nu = 2),
                       tau = list(alpha = alpha_tau, nu = 3),
                       gamma = list(mean = 0, var = 1e1),
                       beta = list(mean = 0, var = 1e1),
                       alpha = list(a = 3, b = 1))

# MCMC simulation
system.time(res_BMR <- bayesmr_mix_het(hdlcad_data, control = control, prior = prior))

res_BMR@results[[1]]@accept
plot(res_BMR@results[[1]]@gamma.chain, type = "l")
lines(res_BMR@results[[2]]@gamma.chain, type = "l", col = gray(.6))
lines(res_BMR@results[[3]]@gamma.chain, type = "l", col = gray(.9))
plot(res_BMR@results[[1]]@alpha.chain, type = "l")
lines(res_BMR@results[[2]]@alpha.chain, type = "l", col = gray(.6))
lines(res_BMR@results[[3]]@alpha.chain, type = "l", col = gray(.9))
plot(res_BMR@results[[1]]@psi.chain, type = "l")
lines(res_BMR@results[[2]]@psi.chain, type = "l", col = gray(.6))
lines(res_BMR@results[[3]]@psi.chain, type = "l", col = gray(.9))
plot(res_BMR@results[[1]]@tau.chain, type = "l")
lines(res_BMR@results[[2]]@tau.chain, type = "l", col = gray(.6))
lines(res_BMR@results[[3]]@tau.chain, type = "l", col = gray(.9))
acf(res_BMR@results[[1]]@tau.chain)
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
  
  simil_probs <- simil_mat/(S - burnin)
  
  return(simil_probs)
}

psm1 <- post_sim(res_BMR@results[[1]]@xi.chain, burnin = 0)
psm2 <- post_sim(res_BMR@results[[2]]@xi.chain, burnin = 0)
psm3 <- post_sim(res_BMR@results[[3]]@xi.chain, burnin = 0)
plot(psm1, psm2); abline(0, 1, lty = 2, col = "gray") ## good agreement
plot(psm1, psm3); abline(0, 1, lty = 2, col = "gray") ## good agreement
plot(psm2, psm3); abline(0, 1, lty = 2, col = "gray") ## good agreement

library(mcclust.ext)
library(mclust)
library(fields)

## Point estimate of the partition
## Dahl (2022) JCGS

library(salso)

## Variation of Information and Binder Loss

vi <- salso(x = res_BMR@results[[1]]@xi.chain, loss = "VI")
table(vi)
bl <- salso(x = res_BMR@results[[1]]@xi.chain, loss = "binder")
table(bl)
touse <- bl

## Wade Gaharamani BA 2018

# vi =  minVI(psm1)$cl
# table(vi)

ordered_ind <- unlist(sapply(length(table(touse)):1,
                                function(k) c(which(touse == k))))

simil_ordered <- psm1[ordered_ind, ordered_ind]

colori <- colorRampPalette(c('white','black'))

p <- ncol(simil_ordered)

par(mar = c(5,5,2,2), mgp = c(3, 1, 0))
image.plot(simil_ordered, col = colori(100),
           zlim = c(0,1), cex.sub = 3, axes = FALSE,
           horizontal = FALSE, legend.shrink = 1,
           xlab = "SNP (j)", ylab = "SNP (j')",
           cex.axis = 4 )
axis(1, at = seq(0, 1, length = p), lab = 1:p, las = 2)
axis(2, at = seq(0, 1, length = p), lab = 1:p, las = 2)

plot(sapply(1:nrow(res_BMR@results[[1]]@xi.chain),
            function(s) length(unique(res_BMR@results[[1]]@xi.chain[s, ]))),
     ylab = "Number of clusters", type = "p", xlab = "MCMC iteration")

S <- nrow(res_BMR@results[[1]]@beta.chain)

# Recover SNP-specific causal effects

beta_snps <- sapply(1:S, function(s) res_BMR@results[[1]]@beta.chain[s, ][res_BMR@results[[1]]@xi.chain[s, ]])
boxplot(t(beta_snps[, 1:S]), outline = FALSE)

# arranged according to the estimated K_start

boxplot(t(beta_snps[ordered_ind, 1:S]), outline = FALSE,
        col = c(rep("lightblue", sum(touse == 1)), rep("pink", sum(touse == 2))))

# Sorted in a more meaningful way

means <- colMeans(t(beta_snps[, 1:S]))
ord <- order(means)

hdl_cad$SNP[ord]
plot(means[ord])

par(mar = c(8, 8, 2, 2), mgp = c(6, 1, 0))
image.plot(psm1[ord, ord], col = colori(100), zlim = c(0,1), cex.sub = 3, 
           axes = FALSE, horizontal = FALSE, legend.shrink = 1,
           xlab = "SNP (j)", ylab = "SNP (j')", cex.axis = 4 )
axis(1, at = seq(0, 1, length = p), lab = hdl_cad$SNP[ord], las = 2)
axis(2, at = seq(0, 1, length = p), lab = hdl_cad$SNP[ord], las = 2)

###

cols <- ifelse(touse == 1, "lightblue", "pink")
beta_post_ord <- t(beta_snps[ord, 1:S])
colnames(beta_post_ord) <- hdl_cad$SNP[ord]
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
