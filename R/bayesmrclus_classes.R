### bayesmr_data ###

#' Input data for a BayesMR analysis.
#'
#' Stores SNP-level summary statistics for a two-sample summary-data Mendelian
#' randomization analysis.
#'
#' @slot data A data frame with four columns: `beta_exposure`, `beta_outcome`,
#'   `se_exposure`, and `se_outcome`. Row names are SNP identifiers taken from
#'   the input `SNP` column. When `reorientation = TRUE` in the initializer,
#'   rows with negative `beta_exposure` are reoriented by changing the sign of
#'   both `beta_exposure` and `beta_outcome`.
#' @slot n A length-one numeric value giving the number of SNPs.
#'
#' @name bayesmr_data-class
#' @rdname bayesmr_data-class
#' @aliases bayesmr_data
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @examples
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' show(dat)
#' summary(dat)
#'
#' @exportClass bayesmr_data
setClass(Class = "bayesmr_data",
  slots = c(
    data = "data.frame",
    n = "numeric"
  )
)

#' Initialize a `bayesmr_data` object.
#'
#' @param .Object Prototype object from class [`bayesmr_data-class`].
#' @param data A data frame containing the columns `SNP`, `beta_exposure`,
#'   `beta_outcome`, `se_exposure`, and `se_outcome`. Dot variants of the four
#'   numeric column names, namely `beta.exposure`, `beta.outcome`, `se.exposure`,
#'   and `se.outcome`, are also accepted and renamed internally.
#' @param reorientation A length-one logical value. If `TRUE`, rows with negative
#'   `beta_exposure` are reoriented by changing the sign of both `beta_exposure`
#'   and `beta_outcome`.
#'
#' @return An initialized [`bayesmr_data-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_data-method
#' @aliases bayesmr_data-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "bayesmr_data",
  function(
    .Object,
    data = data.frame(),
    reorientation = TRUE
  )
  {
    if (any(duplicated(data[, "SNP"])))
      stop("the data set contains duplicated SNPs.")
    cols <- c("beta_exposure", "beta_outcome", "se_exposure", "se_outcome")
    cols2 <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
    if (!all(cols %in% colnames(data))) {
      colnames(data)[match(cols2, colnames(data))] <- cols
    }
    data_tmp <- data[, cols]
    if (any(match(colnames(data), "SNP"), na.rm = TRUE))
      rownames(data_tmp) <- data[, "SNP"]
    if (reorientation) {
      flip <- data_tmp$beta_exposure < 0
      data_tmp$beta_exposure[flip] <- -data_tmp$beta_exposure[flip]
      data_tmp$beta_outcome[flip]  <- -data_tmp$beta_outcome[flip]
    }
    .Object@data <- data_tmp
    .Object@n <- nrow(data_tmp)
    .Object
  }
)

#' Show a `bayesmr_data` object.
#'
#' @param object An object of class [`bayesmr_data-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_data-method
#' @aliases bayesmr_data-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_data",
  function(object) {
    cat("Observed data to use in a BayesMRclus analysis\n")
    cat("Number of SNPs (n):", object@n, "\n")
  }
)

#' Summarize a `bayesmr_data` object.
#'
#' @param object An object of class [`bayesmr_data-class`].
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   summary of the observed data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_data-method
#' @aliases bayesmr_data-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_data",
    function(object, ...) {
      show(object)
      cat("Observed data:\n")
      print(summary(object@data))
      invisible(NULL)
    }
)

#' Plot a `bayesmr_data` object.
#'
#' Produces a scatter plot of SNP-exposure effect estimates against SNP-outcome
#' effect estimates, optionally with horizontal and vertical standard-error bars.
#'
#' @param x An object of class [`bayesmr_data-class`].
#' @param se A length-one logical value indicating whether standard-error bars
#'   should be drawn.
#' @param colors A character vector of colors. This argument is currently accepted
#'   for compatibility but is not used by the plotting code.
#' @param font A numeric value or vector controlling the font for text. This
#'   argument is currently accepted for compatibility but is not used by the
#'   plotting code.
#' @param cex.font A numeric value or vector controlling text expansion. This
#'   argument is currently accepted for compatibility but is not used by the
#'   plotting code.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return Invisibly returns `NULL`; called for its side effect of producing a base
#'   R plot.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_data-method
#' @aliases bayesmr_data-plot
#'
#' @exportMethod plot
#'
#' @examples
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#' plot(dat)
setMethod("plot",
  signature(x = "bayesmr_data"),
  function(x, se = TRUE, colors = c("white", "black"), font = NA, cex.font = NA, ...) {
    D <- x@data
    n <- x@n
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(par(opar))
    graphics::par(mar = c(4, 5, 2, 1) + 0.1, oma = c(1, 0, 0, 0))
    if (se) {
      plot(D[, 1], D[, 2], pch = 20,
        xlab = expression(paste("SNP-exposure effects, ", hat(gamma)[j])),
        ylab = expression(paste("SNP-outcome effects, ", hat(Gamma)[j])),
        xlim = c(min(D[, 1] - D[, 3]), max(D[, 1] + D[, 3])),
        ylim = c(min(D[, 2] - D[, 4]), max(D[, 2] + D[, 4])))
      arrows(D[, 1], D[, 2] - D[, 4], D[, 1], D[, 2] + D[, 4], angle = 0, code = 3, length = 0, lwd = 0.3)
      arrows(D[, 1] - D[, 3], D[, 2], D[, 1] + D[, 3], D[, 2], angle = 0, code = 3, length = 0, lwd = 0.3)
    }
    else {
      plot(D[, 1], D[, 2], pch = 20,
        xlab = expression(paste("SNP-exposure effects, ", hat(gamma)[j])),
        ylab = expression(paste("SNP-outcome effects, ", hat(Gamma)[j])))
      }
    abline(h = 0, lty = 2, col = "gray")
    abline(v = 0, lty = 2, col = "gray")
    points(D[, 1], D[, 2], pch = 20)
  }
)

# ### bayesmr_fit ###

#' Fitted BayesMR model with fixed heterogeneity.
#'
#' Stores the output from a single MCMC chain for a BayesMR model without mixture
#' clustering and with fixed heterogeneity.
#'
#' @slot gamma.chain An array containing posterior draws of the exposure effect
#'   parameter `gamma`.
#' @slot beta.chain An array containing posterior draws of the causal effect
#'   parameter `beta`.
#' @slot accept A numeric vector containing Metropolis-Hastings acceptance
#'   information.
#' @slot data A data frame containing the analyzed summary data.
#' @slot dens A list with log-density components, typically `loglik`, `logprior`,
#'   and `logpost`.
#' @slot control A named list of MCMC control settings.
#' @slot prior A named list of prior hyperparameters.
#' @slot dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @name bayesmr_fit-class
#' @rdname bayesmr_fit-class
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_fit
setClass(Class = "bayesmr_fit",
	slots = c(
		gamma.chain = "array",
		beta.chain = "array",
		accept = "numeric",
		data = "data.frame",
		dens = "list",
		control = "list",
    prior = "list",
		dim = "list"
	)
)

#' Initialize a `bayesmr_fit` object.
#'
#' @param .Object Prototype object from class [`bayesmr_fit-class`].
#' @param gamma.chain An array containing posterior draws of `gamma`.
#' @param beta.chain An array containing posterior draws of `beta` or, for mixture models, cluster-specific causal effects.
#' @param accept A numeric vector or matrix containing Metropolis-Hastings acceptance information.
#' @param data A data frame containing the analyzed summary data.
#' @param dens A list with log-density components, typically `loglik`, `logprior`, and `logpost`.
#' @param control A named list of MCMC control settings.
#' @param prior A named list of prior hyperparameters.
#' @param dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @return An initialized [`bayesmr_fit-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_fit-method
#' @aliases bayesmr_fit-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize",
  "bayesmr_fit",
		function(
			.Object,
			gamma.chain = array(),
			beta.chain = array(),
			accept = numeric(),
			data = data.frame(),
			dens = list(),
			control = list(),
      prior = list(),
			dim = list()
		)
		{
			.Object@gamma.chain <- gamma.chain
			.Object@beta.chain <- beta.chain
			.Object@accept <- accept
			.Object@data <- data
			.Object@dens <- dens
      .Object@control <- control
      .Object@prior <- prior
			.Object@dim <- dim
			.Object
		}
)

#' Show a `bayesmr_fit` object.
#'
#' @param object An object of class [`bayesmr_fit-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_fit-method
#' @aliases bayesmr_fit-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_fit",
  function(object) {
    cat("Bayesian Two-Sample Summary Data simulated chain\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_fit` object.
#'
#' @param object An object of class [`bayesmr_fit-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_fit-method
#' @aliases bayesmr_fit-summary
#'
#' @exportMethod summary
setMethod("summary",
	"bayesmr_fit",
    function(object, include.burnin = FALSE, ...) {
      control <- object@control

      n <- object@dim[["n"]]

      res.coda <- bayesmr_fit_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)

      out <- summary(res.coda)

      return(out)
    }
)

#' Subset an object.
#'
#' Generic function for subsetting BayesMR fitted objects by parameter names.
#'
#' @param x An object to subset.
#' @param ... Further arguments passed to methods.
#'
#' @return The return value depends on the method.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @export
setGeneric("subset", function(x, ...) standardGeneric("subset"))

#' Subset MCMC parameters from a `bayesmr_fit` object.
#'
#' @param x An object of class [`bayesmr_fit-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_fit-method
#' @aliases bayesmr_fit-subset
#'
#' @export
setMethod("subset",
  "bayesmr_fit",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc <- bayesmr_fit_to_mcmc(x, include.burnin = TRUE, verbose = FALSE)

    parnames <- colnames(x_mcmc)
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    
    out <- x_mcmc[, pars]
    
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_fit` object.
#'
#' @param x An object of class [`bayesmr_fit-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_fit-method
#' @aliases bayesmr_fit-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_fit"),
  function(x, what = "trace", pars = character(), regex_pars = "lambda", include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
    
    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc <- bayesmr_fit_to_mcmc(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@control

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    p
  }
)

#' List of fitted BayesMR fixed-heterogeneity chains.
#'
#' Stores the output from one or more MCMC chains for a BayesMR model without
#' mixture clustering and with fixed heterogeneity.
#'
#' @slot results A list of [`bayesmr_fit-class`] objects, one for each simulated
#'   chain.
#'
#' @name bayesmr_fit_list-class
#' @rdname bayesmr_fit_list-class
#' @aliases bayesmr_fit_list
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_fit_list
setClass(Class = "bayesmr_fit_list",
  slots = c(
    results = "list"
  )
)

#' Initialize a `bayesmr_fit_list` object.
#'
#' @param .Object Prototype object from class [`bayesmr_fit_list-class`].
#' @param results A list of fitted chain objects.
#'
#' @return An initialized [`bayesmr_fit_list-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_fit_list-method
#' @aliases bayesmr_fit_list-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "bayesmr_fit_list",
  function(
    .Object,
    results = list()
  )
  {
    .Object@results <- results
    .Object
  }
)

#' Show a `bayesmr_fit_list` object.
#'
#' @param object An object of class [`bayesmr_fit_list-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_fit_list-method
#' @aliases bayesmr_fit_list-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_fit_list",
  function(object) {
    cat("List of Bayesian Two-Sample Summary Data simulated chains\n")
    cat("Number of simulated chains:", length(object@results), "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_fit_list` object.
#'
#' @param object An object of class [`bayesmr_fit_list-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_fit_list-method
#' @aliases bayesmr_fit_list-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_fit_list",
    function(object, include.burnin = FALSE, ...) {
      control <- object@results[[1]]@control
      nchains <- control[["nchains"]]

      n <- object@results[[1]]@dim[["n"]]

      res.coda <- bayesmr_fit_list_to_mcmc.list(object, include.burnin = include.burnin, verbose = FALSE)

      out <- summary(res.coda)

      return(out)
    }
)

#' Subset MCMC parameters from a `bayesmr_fit_list` object.
#'
#' @param x An object of class [`bayesmr_fit_list-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc.list` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_fit_list-method
#' @aliases bayesmr_fit_list-subset
#'
#' @export
setMethod("subset",
  "bayesmr_fit_list",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc.list <- bayesmr_fit_list_to_mcmc.list(x, include.burnin = TRUE, verbose = FALSE)

    parnames <- colnames(x_mcmc.list[[1]])
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    
    out <- coda::mcmc.list(lapply(x_mcmc.list, function(chain) chain[, pars]))
    
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_fit_list` object.
#'
#' @param x An object of class [`bayesmr_fit_list-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_fit_list-method
#' @aliases bayesmr_fit_list-plot
#'
#' @exportMethod plot
setMethod("plot",
	signature(x = "bayesmr_fit_list"),
	function(x, what = "trace", pars = character(), regex_pars = "lambda", include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
		
    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc.list <- bayesmr_fit_list_to_mcmc.list(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@results[[1]]@control

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- control[["nchains"]]*floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    p
	}
)

### bayesmr_het_fit ###

#' Fitted BayesMR model with random heterogeneity.
#'
#' Stores the output from a single MCMC chain for a BayesMR model without mixture
#' clustering and with random heterogeneity.
#'
#' @slot gamma.chain An array containing posterior draws of `gamma`.
#' @slot beta.chain An array containing posterior draws of `beta`.
#' @slot psi.chain An array containing posterior draws of `psi`.
#' @slot tau.chain An array containing posterior draws of `tau`.
#' @slot accept A matrix containing Metropolis-Hastings acceptance information.
#' @slot data A data frame containing the analyzed summary data.
#' @slot dens A list with log-density components, typically `loglik`, `logprior`,
#'   and `logpost`.
#' @slot control A named list of MCMC control settings.
#' @slot prior A named list of prior hyperparameters.
#' @slot dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @name bayesmr_het_fit-class
#' @rdname bayesmr_het_fit-class
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_het_fit
setClass(Class = "bayesmr_het_fit",
  slots = c(
    gamma.chain = "array",
    beta.chain = "array",
    psi.chain = "array",
    tau.chain = "array",
    accept = "matrix",
    data = "data.frame",
    dens = "list",
    control = "list",
    prior = "list",
    dim = "list"
  )
)

#' Initialize a `bayesmr_het_fit` object.
#'
#' @param .Object Prototype object from class [`bayesmr_het_fit-class`].
#' @param gamma.chain An array containing posterior draws of `gamma`.
#' @param beta.chain An array containing posterior draws of `beta` or, for mixture models, cluster-specific causal effects.
#' @param psi.chain An array containing posterior draws of `psi`.
#' @param tau.chain An array containing posterior draws of `tau`.
#' @param accept A numeric vector or matrix containing Metropolis-Hastings acceptance information.
#' @param data A data frame containing the analyzed summary data.
#' @param dens A list with log-density components, typically `loglik`, `logprior`, and `logpost`.
#' @param control A named list of MCMC control settings.
#' @param prior A named list of prior hyperparameters.
#' @param dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @return An initialized [`bayesmr_het_fit-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_het_fit-method
#' @aliases bayesmr_het_fit-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize",
  "bayesmr_het_fit",
    function(
      .Object,
      gamma.chain = array(),
      beta.chain = array(),
      psi.chain = array(),
      tau.chain = array(),
      accept = matrix(),
      data = data.frame(),
      dens = list(),
      control = list(),
      prior = list(),
      dim = list()
    )
    {
      .Object@gamma.chain <- gamma.chain
      .Object@beta.chain <- beta.chain
      .Object@psi.chain <- psi.chain
      .Object@tau.chain <- tau.chain
      .Object@accept <- accept
      .Object@data <- data
      .Object@dens <- dens
      .Object@control <- control
      .Object@prior <- prior
      .Object@dim <- dim
      .Object
    }
)

#' Show a `bayesmr_het_fit` object.
#'
#' @param object An object of class [`bayesmr_het_fit-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_het_fit-method
#' @aliases bayesmr_het_fit-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_het_fit",
  function(object) {
    cat("Bayesian Two-Sample Summary Data simulated chain\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_het_fit` object.
#'
#' @param object An object of class [`bayesmr_het_fit-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_het_fit-method
#' @aliases bayesmr_het_fit-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_het_fit",
    function(object, include.burnin = FALSE, ...) {
      control <- object@control

      n <- object@dim[["n"]]

      res.coda <- bayesmr_het_fit_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)

      out <- summary(res.coda)

      return(out)
    }
)

#' Subset MCMC parameters from a `bayesmr_het_fit` object.
#'
#' @param x An object of class [`bayesmr_het_fit-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_het_fit-method
#' @aliases bayesmr_het_fit-subset
#'
#' @export
setMethod("subset",
  "bayesmr_het_fit",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc <- bayesmr_het_fit_to_mcmc(x, include.burnin = TRUE, verbose = FALSE)

    parnames <- colnames(x_mcmc)
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    
    out <- x_mcmc[, pars]
    
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_het_fit` object.
#'
#' @param x An object of class [`bayesmr_het_fit-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_het_fit-method
#' @aliases bayesmr_het_fit-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_het_fit"),
  function(x, what = "trace", pars = character(), regex_pars = "lambda", include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
    
    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc <- bayesmr_het_fit_to_mcmc(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@control

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    p
  }
)

#' List of fitted BayesMR random-heterogeneity chains.
#'
#' Stores the output from one or more MCMC chains for a BayesMR model without
#' mixture clustering and with random heterogeneity.
#'
#' @slot results A list of [`bayesmr_het_fit-class`] objects, one for each
#'   simulated chain.
#'
#' @name bayesmr_het_fit_list-class
#' @rdname bayesmr_het_fit_list-class
#' @aliases bayesmr_het_fit_list
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_het_fit_list
setClass(Class = "bayesmr_het_fit_list",
  slots = c(
    results = "list"
  )
)

#' Initialize a `bayesmr_het_fit_list` object.
#'
#' @param .Object Prototype object from class [`bayesmr_het_fit_list-class`].
#' @param results A list of fitted chain objects.
#'
#' @return An initialized [`bayesmr_het_fit_list-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_het_fit_list-method
#' @aliases bayesmr_het_fit_list-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "bayesmr_het_fit_list",
  function(
    .Object,
    results = list()
  )
  {
    .Object@results <- results
    .Object
  }
)

#' Show a `bayesmr_het_fit_list` object.
#'
#' @param object An object of class [`bayesmr_het_fit_list-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_het_fit_list-method
#' @aliases bayesmr_het_fit_list-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_het_fit_list",
  function(object) {
    cat("List of Bayesian Two-Sample Summary Data simulated chains\n")
    cat("Number of simulated chains:", length(object@results), "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_het_fit_list` object.
#'
#' @param object An object of class [`bayesmr_het_fit_list-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_het_fit_list-method
#' @aliases bayesmr_het_fit_list-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_het_fit_list",
    function(object, include.burnin = FALSE, ...) {
      control <- object@results[[1]]@control
      nchains <- control[["nchains"]]

      n <- object@results[[1]]@dim[["n"]]

      res.coda <- bayesmr_het_fit_list_to_mcmc.list(object, include.burnin = include.burnin, verbose = FALSE)

      out <- summary(res.coda)

      return(out)
    }
)

#' Subset MCMC parameters from a `bayesmr_het_fit_list` object.
#'
#' @param x An object of class [`bayesmr_het_fit_list-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc.list` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_het_fit_list-method
#' @aliases bayesmr_het_fit_list-subset
#'
#' @export
setMethod("subset",
  "bayesmr_het_fit_list",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc.list <- bayesmr_het_fit_list_to_mcmc.list(x, include.burnin = TRUE, verbose = FALSE)

    parnames <- colnames(x_mcmc.list[[1]])
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    
    out <- coda::mcmc.list(lapply(x_mcmc.list, function(chain) chain[, pars]))
    
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_het_fit_list` object.
#'
#' @param x An object of class [`bayesmr_het_fit_list-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_het_fit_list-method
#' @aliases bayesmr_het_fit_list-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_het_fit_list"),
  function(x, what = "trace", pars = character(), regex_pars = "lambda", include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
    
    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc.list <- bayesmr_het_fit_list_to_mcmc.list(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@results[[1]]@control

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- control[["nchains"]]*floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    p
  }
)
























### EVERYTHING BELOW IS STILL TO COMPLETE ###

### bayesmr_mix_fit ###

#' Fitted BayesMR mixture model with fixed heterogeneity.
#'
#' Stores the output from a single MCMC chain for a BayesMR mixture model with
#' cluster-specific causal effects and fixed heterogeneity.
#'
#' @slot gamma.chain An array containing posterior draws of `gamma`.
#' @slot beta.chain An array containing posterior draws of the cluster-specific
#'   causal effects, with one row per stored MCMC iteration and one column per
#'   SNP.
#' @slot xi.chain An array containing posterior draws of the cluster allocation
#'   indicators, with one row per stored MCMC iteration and one column per SNP.
#' @slot alpha.chain An array containing posterior draws of the mixture
#'   concentration parameter `alpha`.
#' @slot accept A numeric vector containing Metropolis-Hastings acceptance
#'   information.
#' @slot data A data frame containing the analyzed summary data.
#' @slot dens A list with log-density components, typically `loglik`, `logprior`,
#'   and `logpost`.
#' @slot control A named list of MCMC control settings.
#' @slot prior A named list of prior hyperparameters.
#' @slot dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @name bayesmr_mix_fit-class
#' @rdname bayesmr_mix_fit-class
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_mix_fit
setClass(Class = "bayesmr_mix_fit",
  slots = c(
    gamma.chain = "array",
    beta.chain = "array",
    xi.chain = "array",
    alpha.chain = "array",
    accept = "numeric",
    data = "data.frame",
    dens = "list",
    control = "list",
    prior = "list",
    dim = "list"
  )
)

#' Initialize a `bayesmr_mix_fit` object.
#'
#' @param .Object Prototype object from class [`bayesmr_mix_fit-class`].
#' @param gamma.chain An array containing posterior draws of `gamma`.
#' @param beta.chain An array containing posterior draws of `beta` or, for mixture models, cluster-specific causal effects.
#' @param xi.chain An array containing posterior draws of cluster allocation indicators.
#' @param alpha.chain An array containing posterior draws of the mixture concentration parameter `alpha`.
#' @param accept A numeric vector or matrix containing Metropolis-Hastings acceptance information.
#' @param data A data frame containing the analyzed summary data.
#' @param dens A list with log-density components, typically `loglik`, `logprior`, and `logpost`.
#' @param control A named list of MCMC control settings.
#' @param prior A named list of prior hyperparameters.
#' @param dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @return An initialized [`bayesmr_mix_fit-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_mix_fit-method
#' @aliases bayesmr_mix_fit-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize",
  "bayesmr_mix_fit",
    function(
      .Object,
      gamma.chain = array(),
      beta.chain = array(),
      xi.chain = array(),
      alpha.chain = array(),
      accept = numeric(),
      data = data.frame(),
      dens = list(),
      control = list(),
      prior = list(),
      dim = list()
    )
    {
      .Object@gamma.chain <- gamma.chain
      .Object@beta.chain <- beta.chain
      .Object@xi.chain <- xi.chain
      .Object@alpha.chain <- alpha.chain
      .Object@accept <- accept
      .Object@data <- data
      .Object@dens <- dens
      .Object@control <- control
      .Object@prior <- prior
      .Object@dim <- dim
      .Object
    }
)

#' List of fitted BayesMR mixture-model chains.
#'
#' Stores the output from one or more MCMC chains for a BayesMR mixture model with
#' cluster-specific causal effects and fixed heterogeneity.
#'
#' @slot results A list of [`bayesmr_mix_fit-class`] objects, one for each
#'   simulated chain.
#'
#' @name bayesmr_mix_fit_list-class
#' @rdname bayesmr_mix_fit_list-class
#' @aliases bayesmr_mix_fit_list
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_mix_fit_list
setClass(Class = "bayesmr_mix_fit_list",
  slots = c(
    results = "list"
  )
)

#' Initialize a `bayesmr_mix_fit_list` object.
#'
#' @param .Object Prototype object from class [`bayesmr_mix_fit_list-class`].
#' @param results A list of fitted chain objects.
#'
#' @return An initialized [`bayesmr_mix_fit_list-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_mix_fit_list-method
#' @aliases bayesmr_mix_fit_list-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "bayesmr_mix_fit_list",
  function(
    .Object,
    results = list()
  )
  {
    .Object@results <- results
    .Object
  }
)

### bayesmr_mix_het_fit ###

#' Fitted BayesMR mixture model with random heterogeneity.
#'
#' Stores the output from a single MCMC chain for a BayesMR mixture model with
#' cluster-specific causal effects and random heterogeneity.
#'
#' @slot gamma.chain An array containing posterior draws of `gamma`.
#' @slot beta.chain An array containing posterior draws of the cluster-specific
#'   causal effects, with one row per stored MCMC iteration and one column per
#'   SNP.
#' @slot xi.chain An array containing posterior draws of the cluster allocation
#'   indicators, with one row per stored MCMC iteration and one column per SNP.
#' @slot alpha.chain An array containing posterior draws of the mixture
#'   concentration parameter `alpha`.
#' @slot psi.chain An array containing posterior draws of `psi`.
#' @slot tau.chain An array containing posterior draws of `tau`.
#' @slot accept A numeric vector containing Metropolis-Hastings acceptance
#'   information.
#' @slot data A data frame containing the analyzed summary data.
#' @slot dens A list with log-density components, typically `loglik`, `logprior`,
#'   and `logpost`.
#' @slot control A named list of MCMC control settings.
#' @slot prior A named list of prior hyperparameters.
#' @slot dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @name bayesmr_mix_het_fit-class
#' @rdname bayesmr_mix_het_fit-class
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_mix_het_fit
setClass(Class = "bayesmr_mix_het_fit",
  slots = c(
    gamma.chain = "array",
    beta.chain = "array",
    xi.chain = "array",
    alpha.chain = "array",
    psi.chain = "array",
    tau.chain = "array",
    accept = "numeric",
    data = "data.frame",
    dens = "list",
    control = "list",
    prior = "list",
    dim = "list"
  )
)

#' Initialize a `bayesmr_mix_het_fit` object.
#'
#' @param .Object Prototype object from class [`bayesmr_mix_het_fit-class`].
#' @param gamma.chain An array containing posterior draws of `gamma`.
#' @param beta.chain An array containing posterior draws of `beta` or, for mixture models, cluster-specific causal effects.
#' @param xi.chain An array containing posterior draws of cluster allocation indicators.
#' @param alpha.chain An array containing posterior draws of the mixture concentration parameter `alpha`.
#' @param psi.chain An array containing posterior draws of `psi`.
#' @param tau.chain An array containing posterior draws of `tau`.
#' @param accept A numeric vector or matrix containing Metropolis-Hastings acceptance information.
#' @param data A data frame containing the analyzed summary data.
#' @param dens A list with log-density components, typically `loglik`, `logprior`, and `logpost`.
#' @param control A named list of MCMC control settings.
#' @param prior A named list of prior hyperparameters.
#' @param dim A list of model dimensions, including `n`, the number of SNPs.
#'
#' @return An initialized [`bayesmr_mix_het_fit-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_mix_het_fit-method
#' @aliases bayesmr_mix_het_fit-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize",
  "bayesmr_mix_het_fit",
    function(
      .Object,
      gamma.chain = array(),
      beta.chain = array(),
      xi.chain = array(),
      alpha.chain = array(),
      psi.chain = array(),
      tau.chain = array(),
      accept = numeric(),
      data = data.frame(),
      dens = list(),
      control = list(),
      prior = list(),
      dim = list()
    )
    {
      .Object@gamma.chain <- gamma.chain
      .Object@beta.chain <- beta.chain
      .Object@xi.chain <- xi.chain
      .Object@alpha.chain <- alpha.chain
      .Object@psi.chain <- psi.chain
      .Object@tau.chain <- tau.chain
      .Object@accept <- accept
      .Object@data <- data
      .Object@dens <- dens
      .Object@control <- control
      .Object@prior <- prior
      .Object@dim <- dim
      .Object
    }
)

#' List of fitted BayesMR mixture random-heterogeneity chains.
#'
#' Stores the output from one or more MCMC chains for a BayesMR mixture model with
#' cluster-specific causal effects and random heterogeneity.
#'
#' @slot results A list of [`bayesmr_mix_het_fit-class`] objects, one for each
#'   simulated chain.
#'
#' @name bayesmr_mix_het_fit_list-class
#' @rdname bayesmr_mix_het_fit_list-class
#' @aliases bayesmr_mix_het_fit_list
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @exportClass bayesmr_mix_het_fit_list
setClass(Class = "bayesmr_mix_het_fit_list",
  slots = c(
    results = "list"
  )
)

#' Initialize a `bayesmr_mix_het_fit_list` object.
#'
#' @param .Object Prototype object from class [`bayesmr_mix_het_fit_list-class`].
#' @param results A list of fitted chain objects.
#'
#' @return An initialized [`bayesmr_mix_het_fit_list-class`] object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_mix_het_fit_list-method
#' @aliases bayesmr_mix_het_fit_list-initialize
#'
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "bayesmr_mix_het_fit_list",
  function(
    .Object,
    results = list()
  )
  {
    .Object@results <- results
    .Object
  }
)

### bayesmr_mix_fit methods ###

#' Show a `bayesmr_mix_fit` object.
#'
#' @param object An object of class [`bayesmr_mix_fit-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_mix_fit-method
#' @aliases bayesmr_mix_fit-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_mix_fit",
  function(object) {
    cat("Bayesian Two-Sample Summary Data mixture model (single chain)\n")
    cat("Number of SNPs:", object@dim[["n"]], "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_mix_fit` object.
#'
#' @param object An object of class [`bayesmr_mix_fit-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_mix_fit-method
#' @aliases bayesmr_mix_fit-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_mix_fit",
  function(object, include.burnin = FALSE, ...) {
    res.coda <- bayesmr_mix_fit_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)
    out <- summary(res.coda)
    return(out)
  }
)

#' Subset MCMC parameters from a `bayesmr_mix_fit` object.
#'
#' @param x An object of class [`bayesmr_mix_fit-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_mix_fit-method
#' @aliases bayesmr_mix_fit-subset
#'
#' @export
setMethod("subset",
  "bayesmr_mix_fit",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc <- bayesmr_mix_fit_to_mcmc(x, include.burnin = TRUE, verbose = FALSE)
    parnames <- colnames(x_mcmc)
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    out <- x_mcmc[, pars]
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_mix_fit` object.
#'
#' @param x An object of class [`bayesmr_mix_fit-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_mix_fit-method
#' @aliases bayesmr_mix_fit-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_mix_fit"),
  function(x, what = "trace", pars = character(), regex_pars = character(),
           include.burnin = FALSE, combo = NULL, ...) {
    stopifnot(is.character(pars), is.character(regex_pars), is.character(what))

    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc <- bayesmr_mix_fit_to_mcmc(x, include.burnin = include.burnin, verbose = FALSE)
    control <- x@control

    if (what %in% acf_plot_list) {
      if (what == "acf") p <- bayesplot::mcmc_acf(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_acf_bar(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% areas_plot_list) {
      if (what == "areas") p <- bayesplot::mcmc_areas(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_areas_ridges(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% dens_plot_list) {
      if (what == "dens") p <- bayesplot::mcmc_dens(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else if (what == "dens_overlay") p <- bayesplot::mcmc_dens_overlay(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_dens_chains(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hex_plot_list) {
      p <- bayesplot::mcmc_hex(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hist_plot_list) {
      if (what == "hist") p <- bayesplot::mcmc_hist(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% intervals_plot_list) {
      p <- bayesplot::mcmc_intervals(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- floor((control[["burnin"]] + control[["nsim"]]) / control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff / nsample
      if (what == "neff") p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      else p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
    } else if (what %in% pairs_plot_list) {
      p <- bayesplot::mcmc_pairs(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% parcoord_plot_list) {
      p <- bayesplot::mcmc_parcoord(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      else if (what == "recover_intervals") p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      else p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
    } else if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      else p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
    } else if (what %in% scatter_plot_list) {
      p <- bayesplot::mcmc_scatter(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% trace_plot_list) {
      if (what == "trace") p <- bayesplot::mcmc_trace(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_trace_highlight(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% violin_plot_list) {
      p <- bayesplot::mcmc_violin(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what == "combo") {
      if (!is.null(combo))
        p <- bayesplot::mcmc_combo(x = x_mcmc, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      else
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
    }

    p
  }
)

### bayesmr_mix_fit_list methods ###

#' Show a `bayesmr_mix_fit_list` object.
#'
#' @param object An object of class [`bayesmr_mix_fit_list-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_mix_fit_list-method
#' @aliases bayesmr_mix_fit_list-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_mix_fit_list",
  function(object) {
    cat("List of Bayesian Two-Sample Summary Data mixture model chains\n")
    cat("Number of simulated chains:", length(object@results), "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_mix_fit_list` object.
#'
#' @param object An object of class [`bayesmr_mix_fit_list-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_mix_fit_list-method
#' @aliases bayesmr_mix_fit_list-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_mix_fit_list",
  function(object, include.burnin = FALSE, ...) {
    res.coda <- bayesmr_mix_fit_list_to_mcmc.list(object, include.burnin = include.burnin, verbose = FALSE)
    out <- summary(res.coda)
    return(out)
  }
)

#' Subset MCMC parameters from a `bayesmr_mix_fit_list` object.
#'
#' @param x An object of class [`bayesmr_mix_fit_list-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc.list` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_mix_fit_list-method
#' @aliases bayesmr_mix_fit_list-subset
#'
#' @export
setMethod("subset",
  "bayesmr_mix_fit_list",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc.list <- bayesmr_mix_fit_list_to_mcmc.list(x, include.burnin = TRUE, verbose = FALSE)
    parnames <- colnames(x_mcmc.list[[1]])
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    out <- coda::mcmc.list(lapply(x_mcmc.list, function(chain) chain[, pars]))
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_mix_fit_list` object.
#'
#' @param x An object of class [`bayesmr_mix_fit_list-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_mix_fit_list-method
#' @aliases bayesmr_mix_fit_list-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_mix_fit_list"),
  function(x, what = "trace", pars = character(), regex_pars = character(),
           include.burnin = FALSE, combo = NULL, ...) {
    stopifnot(is.character(pars), is.character(regex_pars), is.character(what))

    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc.list <- bayesmr_mix_fit_list_to_mcmc.list(x, include.burnin = include.burnin, verbose = FALSE)
    control <- x@results[[1]]@control

    if (what %in% acf_plot_list) {
      if (what == "acf") p <- bayesplot::mcmc_acf(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_acf_bar(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% areas_plot_list) {
      if (what == "areas") p <- bayesplot::mcmc_areas(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_areas_ridges(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% dens_plot_list) {
      if (what == "dens") p <- bayesplot::mcmc_dens(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else if (what == "dens_overlay") p <- bayesplot::mcmc_dens_overlay(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_dens_chains(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hex_plot_list) {
      p <- bayesplot::mcmc_hex(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hist_plot_list) {
      if (what == "hist") p <- bayesplot::mcmc_hist(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% intervals_plot_list) {
      p <- bayesplot::mcmc_intervals(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- control[["nchains"]] * floor((control[["burnin"]] + control[["nsim"]]) / control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff / nsample
      if (what == "neff") p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      else p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
    } else if (what %in% pairs_plot_list) {
      p <- bayesplot::mcmc_pairs(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% parcoord_plot_list) {
      p <- bayesplot::mcmc_parcoord(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      else if (what == "recover_intervals") p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      else p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
    } else if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      else p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
    } else if (what %in% scatter_plot_list) {
      p <- bayesplot::mcmc_scatter(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% trace_plot_list) {
      if (what == "trace") p <- bayesplot::mcmc_trace(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_trace_highlight(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% violin_plot_list) {
      p <- bayesplot::mcmc_violin(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what == "combo") {
      if (!is.null(combo))
        p <- bayesplot::mcmc_combo(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      else
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
    }

    p
  }
)

### bayesmr_mix_het_fit methods ###

#' Show a `bayesmr_mix_het_fit` object.
#'
#' @param object An object of class [`bayesmr_mix_het_fit-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_mix_het_fit-method
#' @aliases bayesmr_mix_het_fit-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_mix_het_fit",
  function(object) {
    cat("Bayesian Two-Sample Summary Data mixture model with heterogeneity (single chain)\n")
    cat("Number of SNPs:", object@dim[["n"]], "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_mix_het_fit` object.
#'
#' @param object An object of class [`bayesmr_mix_het_fit-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_mix_het_fit-method
#' @aliases bayesmr_mix_het_fit-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_mix_het_fit",
  function(object, include.burnin = FALSE, ...) {
    res.coda <- bayesmr_mix_het_fit_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)
    out <- summary(res.coda)
    return(out)
  }
)

#' Subset MCMC parameters from a `bayesmr_mix_het_fit` object.
#'
#' @param x An object of class [`bayesmr_mix_het_fit-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_mix_het_fit-method
#' @aliases bayesmr_mix_het_fit-subset
#'
#' @export
setMethod("subset",
  "bayesmr_mix_het_fit",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc <- bayesmr_mix_het_fit_to_mcmc(x, include.burnin = TRUE, verbose = FALSE)
    parnames <- colnames(x_mcmc)
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    out <- x_mcmc[, pars]
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_mix_het_fit` object.
#'
#' @param x An object of class [`bayesmr_mix_het_fit-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_mix_het_fit-method
#' @aliases bayesmr_mix_het_fit-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_mix_het_fit"),
  function(x, what = "trace", pars = character(), regex_pars = character(),
           include.burnin = FALSE, combo = NULL, ...) {
    stopifnot(is.character(pars), is.character(regex_pars), is.character(what))

    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc <- bayesmr_mix_het_fit_to_mcmc(x, include.burnin = include.burnin, verbose = FALSE)
    control <- x@control

    if (what %in% acf_plot_list) {
      if (what == "acf") p <- bayesplot::mcmc_acf(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_acf_bar(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% areas_plot_list) {
      if (what == "areas") p <- bayesplot::mcmc_areas(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_areas_ridges(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% dens_plot_list) {
      if (what == "dens") p <- bayesplot::mcmc_dens(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else if (what == "dens_overlay") p <- bayesplot::mcmc_dens_overlay(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_dens_chains(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hex_plot_list) {
      p <- bayesplot::mcmc_hex(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hist_plot_list) {
      if (what == "hist") p <- bayesplot::mcmc_hist(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% intervals_plot_list) {
      p <- bayesplot::mcmc_intervals(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- floor((control[["burnin"]] + control[["nsim"]]) / control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff / nsample
      if (what == "neff") p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      else p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
    } else if (what %in% pairs_plot_list) {
      p <- bayesplot::mcmc_pairs(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% parcoord_plot_list) {
      p <- bayesplot::mcmc_parcoord(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      else if (what == "recover_intervals") p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      else p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
    } else if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      else p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
    } else if (what %in% scatter_plot_list) {
      p <- bayesplot::mcmc_scatter(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% trace_plot_list) {
      if (what == "trace") p <- bayesplot::mcmc_trace(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_trace_highlight(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% violin_plot_list) {
      p <- bayesplot::mcmc_violin(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
    } else if (what == "combo") {
      if (!is.null(combo))
        p <- bayesplot::mcmc_combo(x = x_mcmc, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      else
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
    }

    p
  }
)

### bayesmr_mix_het_fit_list methods ###

#' Show a `bayesmr_mix_het_fit_list` object.
#'
#' @param object An object of class [`bayesmr_mix_het_fit_list-class`].
#'
#' @return Invisibly returns `NULL`; called for its side effect of printing a
#'   short description of the object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_mix_het_fit_list-method
#' @aliases bayesmr_mix_het_fit_list-show
#'
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_mix_het_fit_list",
  function(object) {
    cat("List of Bayesian Two-Sample Summary Data mixture model with heterogeneity chains\n")
    cat("Number of simulated chains:", length(object@results), "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Summarize a `bayesmr_mix_het_fit_list` object.
#'
#' @param object An object of class [`bayesmr_mix_het_fit_list-class`].
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included when converting fitted objects to `coda`
#'   objects before summarizing.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return For fitted objects and fitted-object lists, the result of
#'   `summary.mcmc` or `summary.mcmc.list` via the
#'   corresponding `coda` method. For [`bayesmr_data-class`], invisibly returns
#'   `NULL` after printing a summary of the data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_mix_het_fit_list-method
#' @aliases bayesmr_mix_het_fit_list-summary
#'
#' @exportMethod summary
setMethod("summary",
  "bayesmr_mix_het_fit_list",
  function(object, include.burnin = FALSE, ...) {
    res.coda <- bayesmr_mix_het_fit_list_to_mcmc.list(object, include.burnin = include.burnin, verbose = FALSE)
    out <- summary(res.coda)
    return(out)
  }
)

#' Subset MCMC parameters from a `bayesmr_mix_het_fit_list` object.
#'
#' @param x An object of class [`bayesmr_mix_het_fit_list-class`].
#' @param pars An optional character vector of exact parameter names to keep.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param ... Further arguments passed to or ignored by the method.
#'
#' @return A `coda::mcmc.list` object containing the selected parameters.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases subset,bayesmr_mix_het_fit_list-method
#' @aliases bayesmr_mix_het_fit_list-subset
#'
#' @export
setMethod("subset",
  "bayesmr_mix_het_fit_list",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc.list <- bayesmr_mix_het_fit_list_to_mcmc.list(x, include.burnin = TRUE, verbose = FALSE)
    parnames <- colnames(x_mcmc.list[[1]])
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    out <- coda::mcmc.list(lapply(x_mcmc.list, function(chain) chain[, pars]))
    return(out)
  }
)

#' Plot MCMC diagnostics and summaries for a `bayesmr_mix_het_fit_list` object.
#'
#' @param x An object of class [`bayesmr_mix_het_fit_list-class`].
#' @param what A length-one character vector specifying the plot type. It must be
#'   one of the plot names supported by the package's internal `all_plots_list`
#'   object, such as `"trace"`, `"dens"`, `"hist"`, `"acf"`, `"areas"`,
#'   `"intervals"`, `"rhat"`, `"neff"`, or `"combo"`.
#' @param pars An optional character vector of exact parameter names to plot.
#' @param regex_pars An optional character vector of regular expressions used to
#'   select parameter names.
#' @param include.burnin A length-one logical value indicating whether burn-in
#'   iterations should be included in the plotted draws.
#' @param combo A character vector of plot types passed to
#'   [bayesplot::mcmc_combo()] when `what = "combo"`.
#' @param ... Further arguments passed to the selected `bayesplot` plotting
#'   function.
#'
#' @return A `ggplot` object produced by `bayesplot`.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_mix_het_fit_list-method
#' @aliases bayesmr_mix_het_fit_list-plot
#'
#' @exportMethod plot
setMethod("plot",
  signature(x = "bayesmr_mix_het_fit_list"),
  function(x, what = "trace", pars = character(), regex_pars = character(),
           include.burnin = FALSE, combo = NULL, ...) {
    stopifnot(is.character(pars), is.character(regex_pars), is.character(what))

    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc.list <- bayesmr_mix_het_fit_list_to_mcmc.list(x, include.burnin = include.burnin, verbose = FALSE)
    control <- x@results[[1]]@control

    if (what %in% acf_plot_list) {
      if (what == "acf") p <- bayesplot::mcmc_acf(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_acf_bar(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% areas_plot_list) {
      if (what == "areas") p <- bayesplot::mcmc_areas(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_areas_ridges(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% dens_plot_list) {
      if (what == "dens") p <- bayesplot::mcmc_dens(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else if (what == "dens_overlay") p <- bayesplot::mcmc_dens_overlay(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_dens_chains(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hex_plot_list) {
      p <- bayesplot::mcmc_hex(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% hist_plot_list) {
      if (what == "hist") p <- bayesplot::mcmc_hist(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% intervals_plot_list) {
      p <- bayesplot::mcmc_intervals(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- control[["nchains"]] * floor((control[["burnin"]] + control[["nsim"]]) / control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff / nsample
      if (what == "neff") p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      else p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
    } else if (what %in% pairs_plot_list) {
      p <- bayesplot::mcmc_pairs(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% parcoord_plot_list) {
      p <- bayesplot::mcmc_parcoord(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      else if (what == "recover_intervals") p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      else p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
    } else if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      else p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
    } else if (what %in% scatter_plot_list) {
      p <- bayesplot::mcmc_scatter(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% trace_plot_list) {
      if (what == "trace") p <- bayesplot::mcmc_trace(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      else p <- bayesplot::mcmc_trace_highlight(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what %in% violin_plot_list) {
      p <- bayesplot::mcmc_violin(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
    } else if (what == "combo") {
      if (!is.null(combo))
        p <- bayesplot::mcmc_combo(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      else
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
    }

    p
  }
)
