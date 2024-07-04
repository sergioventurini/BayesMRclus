#' An S4 class to represent the data to use in a BayesMR model.
#'
#' @slot data A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @slot n A length-one character vector representing the number of observations.
#'
#' @name bayesmr_data-class
#' @rdname bayesmr_data-class
#' @aliases bayesmr_data
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2024), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#'
#' @examples
#' showClass("bayesmr_data")
#'
#' @exportClass bayesmr_data
setClass(Class = "bayesmr_data",
  slots = c(
    data = "list",
    n = "numeric"
  )
)

#' Create an instance of the \code{bayesmr_data} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{bayesmr_data}}.
#' @param data A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param n A length-one character vector representing the number of objects
#'   compared by each subject.
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
    data = list(),
    n = numeric()
  )
  {
    .Object@data <- data
    .Object@n <- n
    .Object
  }
)

#' Show an instance of the \code{bayesmr_data} class.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @param object An object of class \code{\link{bayesmr_data}}.
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

#' Provide a summary of a \code{bayesmr_data} class instance.
#'
#' @param object An object of class \code{\link{bayesmr_data}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases summary,bayesmr_data-method
#' @aliases bayesmr_data-summary
#' 
#' @exportMethod summary
setMethod("summary",
  "bayesmr_data",
    function(object) {
      show(object)
      cat("Observed data:\n")
      summary(object@data)
    }
)

#' Provide a graphical summary of a \code{bayesmr_data} class instance.
#'
#' @param x An object of class \code{\link{bayesmr_data}}.
#' @param colors A character vector providing the colors to use in the plot.
#' @param font A length-one numeric vector for the font to use for text.
#'   Can be a vector. \code{NA} values (the default) mean use \code{par("font")}.
#' @param cex.font A length-one numeric vector for the character expansion
#'   factor. \code{NULL} and \code{NA} are equivalent to \code{1.0}. This is an
#'   absolute measure, not scaled by \code{par("cex")} or by setting
#''   \code{par("mfrow")} or \code{par("mfcol")}. Can be a vector.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases plot,bayesmr_data-method
#' @aliases bayesmr_data-plot
#' 
#' @exportMethod plot
#' 
#' @examples
#' data(simdiss)
#' library(bayesplot)
#' cols <- color_scheme_set("brightblue")
#' plot(simdiss, colors = unlist(cols)[c(1, 6)], font = 1, cex.font = 0.75)
setMethod("plot",
  signature(x = "bayesmr_data"),
  function(x, colors = c("white", "black"), font = NA, cex.font = NA, ...) {
    D <- x@data
    n <- x@n
    opar <- graphics::par(no.readonly = TRUE)
    on.exit(par(opar))
    graphics::par(mar = c(4, 5, 2, 1) + 0.1, oma = c(1, 0, 0, 0))
    plot(D[, 1], D[, 2], pch = 19,
      xlab = expression(paste("SNP-exposure effects, ", hat(gamma)[j])),
      ylab = expression(paste("SNP-outcome effects, ", hat(Gamma)[j])))
    abline(h = 0, lty = 2, col = "gray")
    abline(v = 0, lty = 2, col = "gray")
    points(D[, 1], D[, 2], pch = 19)
  }
)

#' An S4 class to represent a BayesMR model.
#'
#' @slot p A length-one character vector representing the number of dimensions
#'   of the latent space to use in the MDS analysis.
#' @slot G A length-one numeric vector representing the number of clusters to
#'   partition the subjects into.
#'
#' @name bayesmr_model-class
#' @rdname bayesmr_model-class
#' @aliases bayesmr_model
#'
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2024), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#'
#' @examples
#' showClass("bayesmr_model")
#'
#' @exportClass bayesmr_model
setClass(Class = "bayesmr_model",
  slots = c(
    p = "numeric",
    G = "numeric"
  )
)

#' Create an instance of the \code{bayesmr_model} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{bayesmr_model}}.
#' @param p A length-one character vector representing the number of dimensions
#'   of the latent space to use in the MDS analysis.
#' @param G A length-one numeric vector representing the number of clusters to
#'   partition the subjects into.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases initialize,bayesmr_model-method
#' @aliases bayesmr_model-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "bayesmr_model",
  function(
    .Object,
    p = numeric(),
    G = numeric()
  )
  {
    .Object@p <- p
    .Object@G <- G
    .Object
  }
)

#' Show an instance of the \code{bayesmr_model} class.
#'
#' @param object An object of class \code{\link{bayesmr_model}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @aliases show,bayesmr_model-method
#' @aliases bayesmr_model-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "bayesmr_model",
  function(object) {
    cat("Bayesian Two-Sample Summary Data Analysis\n")
    cat("Number of latent dimensions (p):", object@p, "\n")
    cat("Number of clusters (G):", object@G, "\n")
  }
)

#' An S4 class to represent the results of fitting BayesMR model.
#'
#' @description
#'   An S4 class to represent the results of fitting BayesMR model using a single
#'   Markov Chain Monte Carlo chain.
#'
#' @slot z.chain An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (untransformed) latent configuration \eqn{Z}.
#' @slot z.chain.p An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (Procrustes-transformed) latent configuration
#'   \eqn{Z}.
#' @slot alpha.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\alpha} parameters.
#' @slot eta.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\eta} parameters.
#' @slot sigma2.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\sigma^2} parameters.
#' @slot lambda.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\lambda} parameters.
#' @slot prob.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership probabilities.
#' @slot x.ind.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership indicators.
#' @slot x.chain An object of class \code{matrix}; posterior draws from
#'   the MCMC algorithm for the cluster membership labels.
#' @slot accept An object of class \code{matrix}; final acceptance rates
#'   for the MCMC algorithm.
#' @slot data An object of class \code{list}; list of observed
#'   dissimilarity matrices.
#' @slot dens An object of class \code{list}; list of log-likelihood,
#'   log-prior and log-posterior values at each iteration of the MCMC simulation.
#' @slot control An object of class \code{list}; list of the control
#'   parameters (number of burnin and sample iterations, number of MCMC chains,
#'   etc.). See \code{\link{bayesmr_control}()} for more information.
#' @slot prior An object of class \code{list}; list of the prior
#'   hyperparameters. See \code{\link{bayesmr_prior}()} for more information.
#' @slot dim An object of class \code{list}; list of dimensions for
#'   the estimated model, i.e. number of objects (\emph{n}), number of latent
#'   dimensions (\emph{p}), number of clusters (\emph{G}), and number of
#'   subjects (\emph{S}).
#' @slot model An object of class \code{\link{bayesmr_model}}.
#'
#' @name bayesmr_fit-class
#' @rdname bayesmr_fit-class
#'
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2024), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#'
#' @examples
#' showClass("bayesmr_fit")
#'
#' @exportClass bayesmr_fit
setClass(Class = "bayesmr_fit",
	slots = c(
		gamma.chain = "array",
		beta.chain = "array",
		accept = "matrix",
		data = "list",
		dens = "list",
		control = "list",
    prior = "list",
		dim = "list",
    model = "bayesmr_model"
	)
)

#' Create an instance of the \code{bayesmr_fit} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{bayesmr_fit}}.
#' @param z.chain An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (untransformed) latent configuration \eqn{Z}.
#' @param z.chain.p An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (Procrustes-transformed) latent configuration
#'   \eqn{Z}.
#' @param alpha.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\alpha} parameters.
#' @param eta.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\eta} parameters.
#' @param sigma2.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\sigma^2} parameters.
#' @param lambda.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\lambda} parameters.
#' @param prob.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership probabilities.
#' @param x.ind.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership indicators.
#' @param x.chain An object of class \code{matrix}; posterior draws from
#'   the MCMC algorithm for the cluster membership labels.
#' @param accept An object of class \code{matrix}; final acceptance rates
#'   for the MCMC algorithm.
#' @param data An object of class \code{list}; list of observed
#'   dissimilarity matrices.
#' @param dens An object of class \code{list}; list of log-likelihood,
#'   log-prior and log-posterior values at each iteration of the MCMC simulation.
#' @param control An object of class \code{list}; list of the control
#'   parameters (number of burnin and sample iterations, number of MCMC chains,
#'   etc.). See \code{\link{bayesmr_control}()} for more information.
#' @param prior An object of class \code{list}; list of the prior
#'   hyperparameters. See \code{\link{bayesmr_prior}()} for more information.
#' @param dim An object of class \code{list}; list of dimensions for
#'   the estimated model, i.e. number of objects (\emph{n}), number of latent
#'   dimensions (\emph{p}), number of clusters (\emph{G}), and number of
#'   subjects (\emph{S}).
#' @param model An object of class \code{\link{bayesmr_model}}.
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
			accept = matrix(),
			data = list(),
			dens = list(),
			control = list(),
      prior = list(),
			dim = list(),
      model = NA
		)
		{
			.Object@gamma.chain <- z.chain
			.Object@beta.chain <- alpha.chain
			.Object@accept <- accept
			.Object@data <- data
			.Object@dens <- dens
      .Object@control <- control
      .Object@prior <- prior
			.Object@dim <- dim
      .Object@model <- model
			.Object
		}
)

#' Show an instance of the \code{bayesmr_fit} class.
#'
#' @param object An object of class \code{\link{bayesmr_fit}}.
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
    cat("Number of latent dimensions (p):", object@model@p, "\n")
    cat("Number of clusters (G):", object@model@G, "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Provide a summary of a \code{bayesmr_fit} class instance.
#'
#' @param object An object of class \code{\link{bayesmr_fit}}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param summary.Z A length-one logical vector. If \code{TRUE} the summary
#'   also includes the latent configuration coordinates.
#' @param ... Further arguments to pass on (currently ignored).
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
      p <- object@dim[["p"]]
      G <- object@dim[["G"]]

      res.coda <- bayesmr_fit_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)

      out <- summary(res.coda)

      return(out)
    }
)

#' Provide a graphical summary of a \code{bayesmr_fit} class instance.
#'
#' @param x An object of class \code{\link{bayesmr_fit}}.
#' @param what A length-one character vector providing the plot type to produce.
#'   Admissible values are those provided by the \pkg{\link{bayesplot}} package,
#'   that is: \code{acf}, \code{areas}, \code{dens}, \code{hex}, \code{hist},
#'   \code{intervals}, \code{neff}, \code{pairs}, \code{parcoord}, \code{recover},
#'   \code{rhat}, \code{scatter}, \code{trace}, \code{violin} or \code{combo}.
#'   In particular, \code{combo} allows to mix different plot types. For more
#'   details see the documentation of the \pkg{\link{bayesplot}} package,
#'   starting from \code{\link[=MCMC-overview]{this overview page}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use for
#'   parameter selection. Can be specified instead of \code{pars} or in addition to
#'   \code{pars}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param combo A character vector providing the plot types to combine (see
#'   \code{\link[bayesplot]{mcmc_combo}}).
#' @param ... Further arguments to pass on.
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

#' An S4 class to represent the results of fitting BayesMR model.
#'
#' @description
#'   An S4 class to represent the results of fitting BayesMR model using multiple
#'   Markov Chain Monte Carlo chains.
#'
#' @slot results An object of class \code{list}; list of \code{bayesmr_fit}
#'   objects corresponding to the parallel MCMC chains simulated during the
#'   estimation.
#'
#' @name bayesmr_fit_list-class
#' @rdname bayesmr_fit_list-class
#' @aliases bayesmr_fit_list
#'
#' @seealso
#' \code{\link{bayesmr_fit}} for more details on the components of each element of
#'   the list.
#'
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2024), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#'
#' @examples
#' showClass("bayesmr_fit_list")
#'
#' @exportClass bayesmr_fit_list
setClass(Class = "bayesmr_fit_list",
  slots = c(
    results = "list"
  )
)

#' Create an instance of the \code{bayesmr_fit_list} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{bayesmr_fit_list}}.
#' @param results A list whose elements are the \code{bayesmr_fit} objects for
#'   each simulated chain.
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

#' Show an instance of the \code{bayesmr_fit_list} class.
#'
#' @param object An object of class \code{\link{bayesmr_fit_list}}.
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
    cat("Number of latent dimensions (p):", object@results[[1]]@model@p, "\n")
    cat("Number of clusters (G):", object@results[[1]]@model@G, "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Provide a summary of a \code{bayesmr_fit_list} class instance.
#'
#' @param object An object of class \code{\link{bayesmr_fit_list}}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param summary.Z A length-one logical vector. If \code{TRUE} the summary
#'   also includes the latent configuration coordinates.
#' @param ... Further arguments to pass on (currently ignored).
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
      p <- object@results[[1]]@dim[["p"]]
      G <- object@results[[1]]@dim[["G"]]

      res.coda <- bayesmr_fit_list_to_mcmc.list(object, include.burnin = include.burnin, verbose = FALSE)

      out <- summary(res.coda)

      return(out)
    }
)

#' Provide a graphical summary of a \code{bayesmr_fit_list} class instance.
#'
#' @param x An object of class \code{\link{bayesmr_fit_list}}.
#' @param what A length-one character vector providing the plot type to produce.
#'   Admissible values are those provided by the \pkg{\link{bayesplot}} package,
#'   that is: \code{acf}, \code{areas}, \code{dens}, \code{hex}, \code{hist},
#'   \code{intervals}, \code{neff}, \code{pairs}, \code{parcoord}, \code{recover},
#'   \code{rhat}, \code{scatter}, \code{trace}, \code{violin} or \code{combo}.
#'   In particular, \code{combo} allows to mix different plot types. For more
#'   details see the documentation of the \pkg{\link{bayesplot}} package,
#'   starting from \code{\link[=MCMC-overview]{this overview page}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use for
#'   parameter selection. Can be specified instead of \code{pars} or in addition to
#'   \code{pars}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param combo A character vector providing the plot types to combine (see
#'   \code{\link[bayesplot]{mcmc_combo}}).
#' @param ... Further arguments to pass on.
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
