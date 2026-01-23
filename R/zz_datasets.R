#' Effect of Body Mass Index (BMI) on Systolic Blood Pressure (SBP)
#'
#' Summary data obtained by combining three genome-wide association studies:
#' \enumerate{
#' \item{BMI-FEM:}{BMI in females by the Genetic Investigation of ANthropometric Traits (GIANT) consortium (sample size: 171977).}
#' \item{BMI-MAL:}{BMI in males in the same study by the GIANT consortium (sample size: 152893)}
#' \item{SBP-UKBB:}{SBP using the United Kingdom BioBank (UKBB) data (sample size: 317754)}
#' }
#'
#' The BMI-FEM dataset is used for SNP selection (column \code{pval.selection}). The BMI-MAL dataset estimates the SNPs' effect on BMI and the SBP-UKBB dataset estimates the SNPs' on SBP.
#'
#' @docType data
#'
#' @usage data(bmi_sbp)
#'
#' @format A \code{data.frame} with 160 rows and 29 variables.
#'
#' @keywords datasets
#'
"bmi_sbp"

#' Effect of HDL choesterol on age-related macular degeneration (AMD)
#'
#' Summary data obtained by combining the following genome-wide association studies:
#' \enumerate{
#' \item{SNP-HDL:}{Kettunen et al. Available from: http://europepmc.org/articles/PMC4814583}
#' \item{SNP-AMD:}{International AMD Genomics Consortium (IAMDGC) Available from: http://europepmc.org/articles/PMC4745342}
#' }
#'
#' The Metabolic Syndrome in Men (METSIM) GWAS has been used for SNP selection.
#'
#' @docType data
#'
#' @usage data(amd_hdl)
#'
#' @format A \code{data.frame} with 27 rows and 4 variables.
#'
#' @keywords datasets
#'
"amd_hdl"

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

#' Effect of Body Mass Index (BMI) on Type-2 Diabetes (T2D)
#'
#' This dataset is created from three genome-wide association studies using the three-sample summary-data MR design
#' \insertCite{zhao2019powerful}{BayesMRclus}:
#' \enumerate{
#' \item \strong{Selection}: \insertCite{akiyama2017genome;textual}{BayesMRclus}
#' \item \strong{Exposure}: \insertCite{locke2015genetic;textual}{BayesMRclus}
#' \item \strong{Outcome}: \insertCite{mahajan2018fine;textual}{BayesMRclus}
#' }
#'
#' The 60 SNPs selected are independent (distance \eqn{\ge 10} mega base pairs, \eqn{R^2 \le 0.001} in a
#' reference panel) and are associated with T2D (p-value less than \eqn{5*10^{-8}}).
#'
#' @references
#' \insertAllCited{}
#'
#' @docType data
#'
#' @usage data(bmi_t2d)
#'
#' @format A \code{data.frame} with 60 rows and 6 variables.
#'
#' @keywords datasets
#'
"bmi_t2d"
