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

#' Educational attainment → Alzheimer’s disease MR dataset
#'
#' Harmonised two-sample Mendelian randomization (MR) dataset
#' constructed using genome-wide significant SNPs for educational
#' attainment as instruments and summary statistics for
#' Alzheimer’s disease as the outcome.
#'
#' Instruments were selected at genome-wide significance
#' (p < 5e-8) and LD-clumped (r2 = 0.001, 10,000 kb window).
#' Proxy SNPs (rsq >= 0.8) were allowed for the outcome GWAS.
#'
#' EXPOSURE: Educational attainment
#'
#' OpenGWAS ID: ieu-a-1239
#'
#' Source study:
#' Lee JJ et al. (2018).
#' Gene discovery and polygenic prediction from a genome-wide
#' association study of educational attainment in 1.1 million individuals.
#' Nature Genetics, 50:1112–1121.
#'
#' GWAS Catalog accession:
#' GCST006250
#'
#' GWAS Catalog URL:
#' https://www.ebi.ac.uk/gwas/studies/GCST006250
#'
#' OpenGWAS dataset page:
#' https://gwas.mrcieu.ac.uk/datasets/ieu-a-1239/
#'
#' OUTCOME: Alzheimer’s disease
#'
#' OpenGWAS ID: ebi-a-GCST90027158
#'
#' Source study:
#' Kunkle BW et al. (2019).
#' Genetic meta-analysis of diagnosed Alzheimer’s disease
#' identifies new risk loci and implicates Aβ, tau,
#' immunity and lipid processing.
#' Nature Genetics, 51:414–430.
#'
#' GWAS Catalog accession:
#' GCST90027158
#'
#' GWAS Catalog URL:
#' https://www.ebi.ac.uk/gwas/studies/GCST90027158
#'
#' OpenGWAS dataset page:
#' https://gwas.mrcieu.ac.uk/datasets/ebi-a-GCST90027158/
#'
#' VARIABLES
#'
#' \describe{
#'   \item{beta_exposure}{Estimated SNP–education association.}
#'   \item{se_exposure}{Standard error of SNP–education association.}
#'   \item{beta_outcome}{Estimated SNP–Alzheimer association (log-odds).}
#'   \item{se_outcome}{Standard error of SNP–Alzheimer association.}
#' }
#'
#' REPRODUCIBILITY
#'
#' The dataset can be regenerated with:
#'
#' \preformatted{
#' educ <- extract_instruments(
#'   outcomes = "ieu-a-1239",
#'   p1 = 5e-08,
#'   clump = TRUE,
#'   r2 = 0.001,
#'   kb = 10000
#' )
#'
#' alzheimer <- extract_outcome_data(
#'   snps = educ$SNP,
#'   outcomes = "ebi-a-GCST90027158",
#'   proxies = TRUE,
#'   rsq = 0.8
#' )
#'
#' dat <- harmonise_data(
#'   exposure_dat = educ,
#'   outcome_dat = alzheimer,
#'   action = 2
#' )
#' }
#'
#' @format A \code{data.frame} with 317 rows and 4 variables.
#'
#' @source OpenGWAS (https://gwas.mrcieu.ac.uk/)
#'
#' @references
#' Lee JJ et al. (2018) Nature Genetics.
#' Kunkle BW et al. (2019) Nature Genetics.
#'
#' @docType data
#'
#' @usage data(educ_alzheimer)
#'
#' @keywords datasets MendelianRandomization GWAS
"educ_alzheimer"
