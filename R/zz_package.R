#' BayesMRclus: Bayesian inference for two-sample summary-data Mendelian
#' randomization under heterogeneity
#'
#' @description
#' The \pkg{BayesMRclus} package provides tools for Bayesian inference in
#' two-sample summary-data Mendelian randomization analyses, allowing for
#' heterogeneity across genetic instruments and, optionally, clustering of
#' instrument-specific causal effects.
#'
#' The main model-fitting functions are:
#' \itemize{
#'   \item [`bayesmr()`]: fixed-heterogeneity model without clustering;
#'   \item [`bayesmr_het()`]: random-heterogeneity model without clustering;
#'   \item [`bayesmr_mix()`]: mixture model with cluster-specific causal effects
#'     and fixed heterogeneity;
#'   \item [`bayesmr_mix_het()`]: mixture model with cluster-specific causal
#'     effects and random heterogeneity.
#' }
#'
#' Supporting functions include [`bayesmr_data-class`] for the input data
#' structure, [`bayesmr_control()`] for MCMC control settings,
#' [`bayesmr_prior()`] for prior hyperparameters, [`bayesmr_init()`] for starting
#' values, and [`run_bnpmr_postprocess()`] for post-processing mixture-model
#' output.
#'
#' @details
#' The package implements MCMC samplers for several BayesMR model variants. The
#' fitted model objects store simulated chains, acceptance information,
#' log-density values, input data, control settings, prior specifications, and
#' model dimensions. Multiple chains can be run serially or in parallel through
#' the options supplied to [`bayesmr_control()`].
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @useDynLib BayesMRclus, .registration = TRUE
#' @import bayesplot
#' @import coda
#' @import ggplot2
#' @import graphics
#' @import grDevices
#' @import mrclust
#' @import MRPATH
#' @import Rdpack
#' @import salso
#' @import stats
#' @import utils
#' @importFrom methods new
#' @importFrom Rcpp evalCpp
#' @importFrom utils prompt .DollarNames
#'
#' @keywords internal
"_PACKAGE"
