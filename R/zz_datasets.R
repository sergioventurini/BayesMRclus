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
