# remotes::install_github('MRCIEU/TwoSampleMR')
# remotes::install_github("mrcieu/ieugwasr")
# remotes::install_github("WSpiller/RadialMR", force = TRUE)
library(TwoSampleMR)
library(ieugwasr)
library(RadialMR)
library(BayesMRclus)

id_instruments <- "ukb-b-19953"  # "ieu-a-2" (smaller --> ~78 SNPs)
id_outcome <- "ieu-a-22"

# check sample overlap
ieugwasr::gwasinfo(id_instruments)$consortium
ieugwasr::gwasinfo(id_outcome)$consortium

# extract data
bmi <- extract_instruments(
  outcomes = id_instruments,
  p1 = 5e-08,
  clump = TRUE,
  r2 = 0.001,
  kb = 10000
)
t2d  <- extract_outcome_data(
  snps = bmi$SNP,
  outcomes = id_outcome,
  proxies = TRUE,
  rsq = 0.8
)
dat_orig <- harmonise_data(
  exposure_dat = bmi,
  outcome_dat = t2d,
  action = 2)
dat_orig <- dat_orig[dat_orig$mr_keep, ]

bmi_t2d <- dat_orig
save(bmi_t2d, file = "/Users/Sergio/dev/BayesMRclus/data/bmi_t2d.rda")

# # checking heterogeneity
# res <- mr(dat_orig)
# mr_scatter_plot(res, dat_orig)

# het <- mr_heterogeneity(dat_orig)
# I2 <- max(0, (het[, "Q"] - (nrow(dat_orig) - 1))/het[, "Q"])

# beta_wald <- dat_orig$beta.outcome/dat_orig$beta.exposure
# hist(beta_wald, breaks = 20)

# radial <- format_radial(dat_orig$beta.exposure,
#                         dat_orig$beta.outcome,
#                         dat_orig$se.exposure,
#                         dat_orig$se.outcome,
#                         dat_orig$SNP)
# ivw.object <- ivw_radial(radial, 0.05, 1, 0.0001, FALSE)
# egg.object <- egger_radial(radial, 0.05, 1, FALSE)
# plot_radial(ivw.object)
# plot_radial(egg.object)

# sd(beta_wald) > median(sqrt((dat_orig$se.outcome^2)/(dat_orig$beta.exposure^2)))

# z_exp <- abs(dat_orig$beta.exposure/dat_orig$se.exposure)
# plot(z_exp, beta_wald, pch = 20,
#      col="darkblue",
#      xlab="Instrument strength |Z_gamma|",
#      ylab="Wald ratio estimate")
# abline(h = 0, lty = 2)
# # add IVW estimate
# ivw <- mr(dat_orig, method_list = "mr_ivw")
# abline(h = ivw$b, col = "red", lwd = 2)
# # smooth trend
# lines(lowess(z_exp, beta_wald), col = "black", lwd = 2)
# precision <- 1/dat_orig$se.outcome
# symbols(z_exp, beta_wald,
#         circles = precision/max(precision),
#         inches = 0.1,
#         fg = "darkblue",
#         add = TRUE)

# # devtools::install_github("cnfoley/mrclust")
# library(ggplot2)
# library(dplyr)
# library(mrclust)

# data("bmi_t2d", package = "BayesMRclus")
# mrclust_res <- mrclust_results(bmi_t2d)
# mrclust_res$p1
# mrclust_res$tab
