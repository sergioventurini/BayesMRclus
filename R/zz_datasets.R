#' Effect of HDL Cholesterol (HDL-C) on Coronary Heart Disease (CHD)
#'
#' This dataset is created from three genome-wide association studies using the three-sample summary-data
#' MR design \insertCite{zhao2019powerful}{BayesMRclus}:
#' \enumerate{
#' \item{Selection}: GWAS of HDL-C by \insertCite{teslovich2010biological;textual}{BayesMRclus}
#' \item{Exposure}: GWAS of lipoprotein subfractions by \insertCite{kettunen2016genome;textual}{BayesMRclus}
#' \item{Outcome}: The CARDIoGRAMplusC4D with 1000 Genome Project imputation GWAS of CAD \insertCite{nikpay2015comprehensive}{BayesMRclus}.
#' }
#'
#' The 151 SNPs selected are independent (distance \eqn{\ge 10} mega base pairs, \eqn{R^2 \le 0.001} in a
#' reference panel) and are associated with at least one plasma lipid trait (the minimum p-value with
#' HDl-C, LDL-C, and triglycerides is less than \eqn{10^{-4}}).
#'
#' @references
#' \insertAllCited{}
#'
#' @docType data
#'
#' @usage data(hdl_chd)
#'
#' @format A \code{data.frame} with 151 rows and 6 variables.
#'
#' @keywords datasets
#'
"hdl_chd"

#' LDL Cholesterol and Coronary Artery Disease GWAS Summary Statistics
#'
#' @description
#' Harmonised two-sample Mendelian randomisation (MR) summary statistics for the
#' causal effect of low-density lipoprotein (LDL) cholesterol on coronary artery
#' disease (CAD), derived from publicly available genome-wide association studies
#' (GWAS) accessed through the OpenGWAS database
#' \insertCite{Elsworth2020}{BayesMRclus}.
#'
#' The dataset was constructed using the \pkg{TwoSampleMR} R package
#' \insertCite{Hemani2018}{BayesMRclus}. Genetic instruments for LDL cholesterol
#' were selected at genome-wide significance and pruned for linkage
#' disequilibrium. Exposure and outcome associations were harmonised with
#' \code{action = 2}; variants passing the \code{mr_keep} harmonisation filter
#' were retained.
#'
#' @format A data frame with 73 rows and 44 variables:
#' \describe{
#'   \item{SNP}{Character. SNP identifier in dbSNP \code{rs} format.}
#'   \item{effect_allele.exposure}{Character. Effect allele in the exposure GWAS.}
#'   \item{other_allele.exposure}{Character. Other allele in the exposure GWAS.}
#'   \item{effect_allele.outcome}{Character. Effect allele in the outcome GWAS.}
#'   \item{other_allele.outcome}{Character. Other allele in the outcome GWAS.}
#'   \item{beta.exposure}{Numeric. SNP association with LDL cholesterol.}
#'   \item{beta.outcome}{Numeric. SNP association with coronary artery disease, expressed as a log-odds ratio.}
#'   \item{eaf.exposure}{Numeric. Effect allele frequency in the exposure GWAS.}
#'   \item{eaf.outcome}{Numeric. Effect allele frequency in the outcome GWAS.}
#'   \item{remove}{Logical. Harmonisation flag indicating variants to remove.}
#'   \item{palindromic}{Logical. Whether the variant is palindromic.}
#'   \item{ambiguous}{Logical. Whether the variant is ambiguous.}
#'   \item{id.outcome}{Character. OpenGWAS outcome identifier.}
#'   \item{chr}{Integer. Chromosome in the outcome GWAS.}
#'   \item{pos}{Integer. Base-pair position in the outcome GWAS.}
#'   \item{se.outcome}{Numeric. Standard error of \code{beta.outcome}.}
#'   \item{samplesize.outcome}{Integer. Outcome GWAS sample size.}
#'   \item{pval.outcome}{Numeric. Outcome GWAS association p-value.}
#'   \item{outcome}{Character. Outcome trait name and identifier.}
#'   \item{originalname.outcome}{Character. Original outcome trait name.}
#'   \item{outcome.deprecated}{Character. Deprecated outcome trait label.}
#'   \item{mr_keep.outcome}{Logical. Outcome-side harmonisation keep flag.}
#'   \item{data_source.outcome}{Character. Outcome data source.}
#'   \item{proxy.outcome}{Logical. Whether the outcome association used a proxy SNP.}
#'   \item{target_snp.outcome}{Character. Target SNP for proxy lookup, if any.}
#'   \item{proxy_snp.outcome}{Character. Proxy SNP used for the outcome, if any.}
#'   \item{target_a1.outcome}{Character. Target effect allele for proxy lookup.}
#'   \item{target_a2.outcome}{Character. Target other allele for proxy lookup.}
#'   \item{proxy_a1.outcome}{Character. Proxy effect allele.}
#'   \item{proxy_a2.outcome}{Character. Proxy other allele.}
#'   \item{id.exposure}{Character. OpenGWAS exposure identifier.}
#'   \item{chr.exposure}{Integer. Chromosome in the exposure GWAS.}
#'   \item{pos.exposure}{Integer. Base-pair position in the exposure GWAS.}
#'   \item{se.exposure}{Numeric. Standard error of \code{beta.exposure}.}
#'   \item{pval.exposure}{Numeric. Exposure GWAS association p-value.}
#'   \item{samplesize.exposure}{Numeric. Exposure GWAS sample size.}
#'   \item{exposure}{Character. Exposure trait name and identifier.}
#'   \item{mr_keep.exposure}{Logical. Exposure-side harmonisation keep flag.}
#'   \item{pval_origin.exposure}{Character. Source of the exposure p-value.}
#'   \item{data_source.exposure}{Character. Exposure data source.}
#'   \item{action}{Integer. Harmonisation action used by \pkg{TwoSampleMR}.}
#'   \item{SNP_index}{Integer. SNP index.}
#'   \item{mr_keep}{Logical. Overall harmonisation keep flag.}
#' }
#'
#' @details
#' \strong{Allele harmonisation.}
#' Exposure and outcome associations are reported with respect to harmonised
#' effect alleles. The exposure associations are not reoriented to be positive;
#' therefore \code{beta.exposure} contains both positive and negative values.
#'
#' @examples
#' data(ldl_cad)
#'
#' wald <- ldl_cad$beta.outcome / ldl_cad$beta.exposure
#' hist(
#'   wald,
#'   breaks = 30,
#'   main = "Wald ratio distribution: LDL cholesterol on CAD",
#'   xlab = "Wald ratio"
#' )
#'
#' w <- 1 / ldl_cad$se.outcome^2
#' beta_ivw <- sum(
#'   w * ldl_cad$beta.outcome * ldl_cad$beta.exposure
#' ) / sum(
#'   w * ldl_cad$beta.exposure^2
#' )
#'
#' cat(sprintf("IVW estimate (log-OR): %.4f\n", beta_ivw))
#'
#' mcmc_dat <- data.frame(
#'   gamma_hat = ldl_cad$beta.exposure,
#'   Gamma_hat = ldl_cad$beta.outcome,
#'   sigma2X   = ldl_cad$se.exposure^2,
#'   sigma2Y   = ldl_cad$se.outcome^2
#' )
"ldl_cad"

#' WHRadjBMI and Fasting Insulin GWAS Summary Statistics
#'
#' @description
#' Harmonised two-sample Mendelian randomisation (MR) summary statistics
#' for the causal effect of waist-to-hip ratio adjusted for body mass
#' index (WHRadjBMI) on fasting blood insulin. Data were derived from
#' publicly available genome-wide association studies (GWAS) accessed
#' through the OpenGWAS database \insertCite{Elsworth2020}{BayesMRclus}
#' and harmonised using the \pkg{TwoSampleMR} R package
#' \insertCite{Hemani2018}{BayesMRclus}.
#'
#' The dataset contains 213 genetic instruments for WHRadjBMI selected
#' at genome-wide significance (\eqn{p < 5 \times 10^{-8}}) and pruned
#' for linkage disequilibrium using \eqn{r^2 < 0.001} within a 10 Mb
#' window. Exposure and outcome associations were harmonised with
#' \code{action = 2}, allowing palindromic variants to be aligned using
#' allele frequency information where possible. Only variants passing
#' the \code{mr_keep} harmonisation filter were retained.
#'
#' This dataset is used as a primary empirical application for evaluating
#' Bayesian nonparametric mixture MR methods in a setting with substantial
#' instrument heterogeneity.
#'
#' @format A data frame with 213 rows and 44 variables:
#' \describe{
#'   \item{SNP}{Character. SNP identifier in dbSNP \code{rs} format.}
#'   \item{effect_allele.exposure}{Character. Effect allele in the exposure GWAS.}
#'   \item{other_allele.exposure}{Character. Other allele in the exposure GWAS.}
#'   \item{effect_allele.outcome}{Character. Effect allele in the outcome GWAS.}
#'   \item{other_allele.outcome}{Character. Other allele in the outcome GWAS.}
#'   \item{beta.exposure}{Numeric. SNP association with WHRadjBMI.}
#'   \item{beta.outcome}{Numeric. SNP association with fasting insulin.}
#'   \item{eaf.exposure}{Numeric. Effect allele frequency in the exposure GWAS.}
#'   \item{eaf.outcome}{Numeric. Effect allele frequency in the outcome GWAS.}
#'   \item{remove}{Logical. Harmonisation flag indicating variants to remove.}
#'   \item{palindromic}{Logical. Whether the variant is palindromic.}
#'   \item{ambiguous}{Logical. Whether the variant is ambiguous.}
#'   \item{id.outcome}{Character. OpenGWAS outcome identifier.}
#'   \item{chr}{Integer. Chromosome in the outcome GWAS.}
#'   \item{pos}{Integer. Base-pair position in the outcome GWAS.}
#'   \item{se.outcome}{Numeric. Standard error of \code{beta.outcome}.}
#'   \item{samplesize.outcome}{Integer. Outcome GWAS sample size.}
#'   \item{pval.outcome}{Numeric. Outcome GWAS association p-value.}
#'   \item{outcome}{Character. Outcome trait name and identifier.}
#'   \item{originalname.outcome}{Character. Original outcome trait name.}
#'   \item{outcome.deprecated}{Character. Deprecated outcome trait label.}
#'   \item{mr_keep.outcome}{Logical. Outcome-side harmonisation keep flag.}
#'   \item{data_source.outcome}{Character. Outcome data source.}
#'   \item{proxy.outcome}{Logical. Whether the outcome association used a proxy SNP.}
#'   \item{target_snp.outcome}{Character. Target SNP for proxy lookup, if any.}
#'   \item{proxy_snp.outcome}{Character. Proxy SNP used for the outcome, if any.}
#'   \item{target_a1.outcome}{Character. Target effect allele for proxy lookup.}
#'   \item{target_a2.outcome}{Character. Target other allele for proxy lookup.}
#'   \item{proxy_a1.outcome}{Character. Proxy effect allele.}
#'   \item{proxy_a2.outcome}{Character. Proxy other allele.}
#'   \item{id.exposure}{Character. OpenGWAS exposure identifier.}
#'   \item{chr.exposure}{Integer. Chromosome in the exposure GWAS.}
#'   \item{pos.exposure}{Integer. Base-pair position in the exposure GWAS.}
#'   \item{se.exposure}{Numeric. Standard error of \code{beta.exposure}.}
#'   \item{pval.exposure}{Numeric. Exposure GWAS association p-value.}
#'   \item{samplesize.exposure}{Integer. Exposure GWAS sample size.}
#'   \item{exposure}{Character. Exposure trait name and identifier.}
#'   \item{mr_keep.exposure}{Logical. Exposure-side harmonisation keep flag.}
#'   \item{pval_origin.exposure}{Character. Source of the exposure p-value.}
#'   \item{data_source.exposure}{Character. Exposure data source.}
#'   \item{action}{Integer. Harmonisation action used by \pkg{TwoSampleMR}.}
#'   \item{SNP_index}{Integer. SNP index.}
#'   \item{mr_keep}{Logical. Overall harmonisation keep flag.}
#' }
#'
#' @details
#' \strong{Exposure GWAS.}
#' WHRadjBMI summary statistics were obtained from the GIANT consortium
#' sex-combined GWAS meta-analysis (\eqn{n = 458{,}349}; OpenGWAS
#' accession \code{ebi-a-GCST90025996})
#' \insertCite{Shungin2015}{BayesMRclus}. WHRadjBMI is a measure of body
#' fat distribution adjusted for overall adiposity.
#'
#' \strong{Outcome GWAS.}
#' Fasting insulin summary statistics were obtained from the MAGIC
#' consortium GWAS (\eqn{n = 51{,}750}; OpenGWAS accession
#' \code{ebi-a-GCST005185}) \insertCite{Manning2012}{BayesMRclus}.
#' The outcome is inverse-normal transformed log-fasting insulin, so
#' effect sizes are not directly interpretable in clinical units.
#'
#' \strong{Instrument selection and harmonisation.}
#' Instruments were extracted using
#' \code{TwoSampleMR::extract_instruments()} with
#' \code{p1 = 5e-8}, \code{clump = TRUE},
#' \code{r2 = 0.001}, and \code{kb = 10000}. Outcome associations were
#' extracted from the fasting insulin GWAS and harmonised using
#' \code{harmonise_data(action = 2)}. Variants failing harmonisation or
#' marked with \code{mr_keep == FALSE} were removed.
#'
#' \strong{Instrument strength.}
#' All retained instruments have strong associations with WHRadjBMI
#' (\eqn{z_j = \hat\gamma_j / s_{\hat\gamma_j} > 5}), corresponding to
#' approximate per-variant \eqn{F}-statistics greater than 25.
#'
#' \strong{Heterogeneity.}
#' This application shows substantial heterogeneity in SNP-specific
#' causal estimates. Preliminary MRClust analysis
#' \insertCite{Foley2021}{BayesMRclus} identified four components:
#' a positive causal cluster, a negative causal cluster, a Null component,
#' and a Junk component. This structure motivates the use of mixture MR
#' models capable of representing heterogeneous causal pathways.
#'
#' \strong{Scientific context.}
#' Genetic instruments for WHRadjBMI capture multiple biological pathways
#' related to adipose storage capacity, visceral fat deposition, insulin
#' resistance, and sex-hormone signalling. Consequently, the relationship
#' between body fat distribution and fasting insulin is expected to show
#' more heterogeneity than benchmark MR examples with a single dominant
#' causal pathway.
#'
#' \strong{Data provenance.}
#' Summary statistics were downloaded from OpenGWAS on April 7, 2026
#' using \pkg{TwoSampleMR} version 0.7.4.
#'
#' @source
#' Exposure GWAS: OpenGWAS accession \code{ebi-a-GCST90025996}.
#' \url{https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90025996/}
#'
#' Outcome GWAS: OpenGWAS accession \code{ebi-a-GCST005185}.
#' \url{https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST005185/}
#'
#' Harmonisation and instrument extraction performed using
#' \pkg{TwoSampleMR}:
#' \url{https://mrcieu.github.io/TwoSampleMR/}
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data(whradjbmi_insulin)
#'
#' # Wald ratios
#' wald <- whradjbmi_insulin$beta.outcome /
#'   whradjbmi_insulin$beta.exposure
#'
#' hist(
#'   wald,
#'   breaks = 30,
#'   main = "Wald ratio distribution: WHRadjBMI on fasting insulin",
#'   xlab = "Wald ratio"
#' )
#'
#' # Inverse-variance weighted estimate
#' w <- 1 / whradjbmi_insulin$se.outcome^2
#' beta_ivw <- sum(
#'   w * whradjbmi_insulin$beta.outcome *
#'     whradjbmi_insulin$beta.exposure
#' ) / sum(
#'   w * whradjbmi_insulin$beta.exposure^2
#' )
#'
#' cat(sprintf("IVW estimate: %.4f\n", beta_ivw))
#'
#' # Prepare input for bayesmr_mixture()
#' mcmc_dat <- data.frame(
#'   gamma_hat = whradjbmi_insulin$beta.exposure,
#'   Gamma_hat = whradjbmi_insulin$beta.outcome,
#'   sigma2X   = whradjbmi_insulin$se.exposure^2,
#'   sigma2Y   = whradjbmi_insulin$se.outcome^2
#' )
"whradjbmi_insulin"
