library(TwoSampleMR)
library(RadialMR)
library(BayesMRclus)

data("DBP_CAD", package = "mrclust")
dbp_cad <- DBP_CAD[, 2:6]
colnames(dbp_cad) <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
save(dbp_cad, file = "/Users/Sergio/dev/BayesMRclus/data/dbp_cad.rda")

# dat_orig <- dbp_cad
# 
# beta_wald <- dat_orig$beta.outcome/dat_orig$beta.exposure
# se_beta_wald <- abs(dat_orig$se.outcome/dat_orig$beta.exposure)
# hist(beta_wald, breaks = 20)
# I2(beta_wald, se_beta_wald)
# 
# radial <- format_radial(dat_orig$beta.exposure,
#                         dat_orig$beta.outcome,
#                         dat_orig$se.exposure,
#                         dat_orig$se.outcome,
#                         dat_orig$SNP)
# ivw.object <- ivw_radial(radial, 0.05, 1, 0.0001, FALSE)
# egg.object <- egger_radial(radial, 0.05, 1, FALSE)
# plot_radial(ivw.object)
# plot_radial(egg.object)
# 
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
# 
# sel <- dat_orig$beta.exposure < 0
# dat_orig$beta.exposure[sel] <- -dat_orig$beta.exposure[sel]
# dat_orig$beta.outcome[sel] <- -dat_orig$beta.outcome[sel]
# 
# mrclust_res <- mrclust_results(dat_orig)
# mrclust_res$p1
# mrclust_res$tab
