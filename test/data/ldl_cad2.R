# remotes::install_github('MRCIEU/TwoSampleMR')
# remotes::install_github("mrcieu/ieugwasr")
# remotes::install_github("WSpiller/RadialMR", force = TRUE)
library(TwoSampleMR)
library(ieugwasr)
library(RadialMR)
library(BayesMRclus)

id_instruments <- "ebi-a-GCST90018961"
id_outcome <- "ieu-a-7"

# check sample overlap
ieugwasr::gwasinfo(id_instruments)$consortium
ieugwasr::gwasinfo(id_outcome)$consortium

# extract data
ldl <- extract_instruments(
  outcomes = id_instruments,
  p1 = 5e-08,
  clump = TRUE,
  r2 = 0.001,
  kb = 10000
)
cad2 <- extract_outcome_data(
  snps = ldl$SNP,
  outcomes = id_outcome,
  proxies = TRUE,
  rsq = 0.8
)
dat_orig <- harmonise_data(
  exposure_dat = ldl,
  outcome_dat = cad2,
  action = 2)

ldl_cad2 <- dat_orig
save(ldl_cad2, file = "/Users/Sergio/dev/BayesMRclus/data/ldl_cad2.rda")

# remove extreme tail outliers
wald <- dat_orig$beta.outcome / dat_orig$beta.exposure
keep <- abs(wald) < 3
dat_orig <- dat_orig[keep, ]

# checking heterogeneity
res <- mr(dat_orig)
mr_scatter_plot(res, dat_orig)

het <- mr_heterogeneity(dat_orig)
I2 <- max(0, (het[, "Q"] - (nrow(dat_orig) - 1))/het[, "Q"])

beta_wald <- dat_orig$beta.outcome/dat_orig$beta.exposure
hist(beta_wald, breaks = 20)

radial <- format_radial(dat_orig$beta.exposure,
                        dat_orig$beta.outcome,
                        dat_orig$se.exposure,
                        dat_orig$se.outcome,
                        dat_orig$SNP)
ivw.object <- ivw_radial(radial, 0.05, 1, 0.0001, FALSE)
egg.object <- egger_radial(radial, 0.05, 1, FALSE)
plot_radial(ivw.object)
plot_radial(egg.object)

sd(beta_wald) > median(sqrt((dat_orig$se.outcome^2)/(dat_orig$beta.exposure^2)))

z_exp <- abs(dat_orig$beta.exposure/dat_orig$se.exposure)
plot(z_exp, beta_wald, pch = 20,
     col="darkblue",
     xlab="Instrument strength |Z_gamma|",
     ylab="Wald ratio estimate")
abline(h = 0, lty = 2)
# add IVW estimate
ivw <- mr(dat_orig, method_list = "mr_ivw")
abline(h = ivw$b, col = "red", lwd = 2)
# smooth trend
lines(lowess(z_exp, beta_wald), col = "black", lwd = 2)
precision <- 1/dat_orig$se.outcome
symbols(z_exp, beta_wald,
        circles = precision/max(precision),
        inches = 0.1,
        fg = "darkblue",
        add = TRUE)

# devtools::install_github("cnfoley/mrclust")
library(ggplot2)
library(dplyr)
library(mrclust)

data("ldl_cad2", package = "BayesMRclus")
mrclust_res <- mrclust_results(ldl_cad2)
mrclust_res$p1
mrclust_res$tab
