#' Estimation of a BayesMR model.
#'
#' \code{bayesmr_mix_het()}, the main function of the package, estimates a BayesMR model
#'   for a given set of \emph{S} dissimilarity matrices.
#'
#' @param data An object of class \code{bayesmr_data} containing the data
#'   to analyze.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution. See
#'   \code{\link{bayesmr_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{bayesmr_prior}()} for more details.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created for the duration of the \code{bayesmr_mix_het()} call.
#' @param post_all A length-one logical vector, which if TRUE applies a further
#'   post-processing to the simulated chains (in case these are more than one).
#' @return A \code{bayesmr_fit_list} object.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bmds}} for Bayesian (metric) multidimensional scaling.
#' @seealso \code{\link{bayesmr_data}} for a description of the data format.
#' @seealso \code{\link{bayesmr_fit_list}} for a description of the elements
#'   included in the returned object.
#' @references
#'   Consonni, G., Venturini, S., Castelletti, F. (2026), "Bayesian Hierarchical Modeling for
#'   Two-Sample Summary-Data Mendelian Randomization under Heterogeneity, working paper.
#' @examples
#' \dontrun{
#' data(simdiss, package = "bayesmr")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 20000
#' nsim <- 10000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random_start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.bayesmr <- bayesmr_mix_het(simdiss, control)
#'
#' summary(sim.bayesmr, include.burnin = FALSE)
#'
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("teal")
#' plot(sim.bayesmr, what = "trace", regex_pars = "eta")
#'
#' z <- bayesmr_get_configuration(sim.bayesmr, chain = 1, est = "mean",
#'   labels = 1:16)
#' summary(z)
#' color_scheme_set("mix-pink-blue")
#' graph <- plot(z, size = 2, size_lbl = 3)
#' graph + panel_bg(fill = "gray90", color = NA)
#' }
#'
#' @importFrom abind abind
#' @export
bayesmr_mix_het <- function(data, control = bayesmr_control(),
  prior = NULL, cl = NULL, post_all = FALSE) {
  
  control <- check_list_na(control, bayesmr_control())
  if (!check_control(control))
    stop("the control list is not correct; see the documentation for more details.")

  nsim <- control[["nsim"]]
  burnin <- control[["burnin"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  threads <- control[["threads"]]
  seed <- control[["seed"]]
  parallel <- control[["parallel"]]
  random_start <- control[["random_start"]]
  if (is.null(control[["K_start"]]) || control[["K_start"]] < 2) {
    control[["K_start"]] <- 2
    message("initial number of clusters set to two.")
  }
  K_start <- control[["K_start"]]
  store.burnin <- control[["store.burnin"]]
  verbose <- control[["verbose"]]

  n <- nrow(data)
  totiter <- burnin + nsim
  
  if (is.null(prior)) {
    prior <- bayesmr_prior()
  } else {
    prior <- check_list_na(prior, bayesmr_prior())
  }
  if (!check_prior(prior)) {
    stop("the prior hyperparameter list is not correct; see the documentation for more details.")
  }

  have_mc <- have_snow <- FALSE
  if (parallel != "no" && threads > 1L) {
    if (parallel == "multicore") {
      have_mc <- (.Platform$OS.type != "windows")
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }
    if (!have_mc && !have_snow) {
      warning("number of cores forced to 1 (i.e. no parallel computing used).")
      threads <- 1L
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }

  use_parallel <- (nchains > 1L) && (have_mc || have_snow)

  ## Preserve caller RNG state completely and restore it on exit
  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  old_rng_kinds <- RNGkind()

  on.exit({
    do.call(RNGkind, as.list(old_rng_kinds))
    if (had_random_seed) {
      assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  ## Derive an internal seed.
  ## - if user supplied one, use it
  ## - otherwise, derive one from the current session RNG so behavior depends on
  ##   the current session state without permanently mutating it outside the function
  internal_seed <- seed
  if (is.null(internal_seed)) {
    internal_seed <- sample.int(.Machine$integer.max - 1L, size = 1L)
  }

  bayesmr_mix_het_fit_parallel <- function(c, data.c, control.c, prior.c, lib) {
    suppressMessages(require(BayesMRclus, lib.loc = lib))
    control.c[["verbose"]] <- FALSE

    start.c <- bayesmr_init(
      data = data.c,
      random_start = control.c[["random_start"]],
      K_start = control.c[["K_start"]]
    )
    
    bayesmr_mix_het_fit(
      data = data.c,
      control = control.c,
      prior = prior.c,
      start = start.c
    )
  }

  bayesmr_mix_het_fit_serial <- function(ch, data.c, control.c, prior.c) {
    if (control.c[["verbose"]] && control.c[["nchains"]] > 1L) {
      message("--- STARTING SIMULATION OF CHAIN ", ch, " OF ", control.c[["nchains"]], " ---")
    }
    
    if (control.c[["verbose"]]) {
      message("Initialization of the algorithm...")
    }
    
    start.c <- bayesmr_init(
      data = data.c,
      random_start = control.c[["random_start"]],
      K_start = control.c[["K_start"]]
    )
    
    fit <- bayesmr_mix_het_fit(
      data = data.c,
      control = control.c,
      prior = prior.c,
      start = start.c
    )
    
    if (control.c[["verbose"]] && control.c[["nchains"]] > 1L) {
      message("--- END OF CHAIN ", ch, " OF ", control.c[["nchains"]], " ---\n")
    }
    
    fit
  }

  # perform MCMC simulation
  if (use_parallel) {
    ## Parallel-safe RNG streams
    RNGkind("L'Ecuyer-CMRG")
    set.seed(internal_seed)

    if (verbose) {
      devout <- ""
      if (.Platform$OS.type != "windows" && !have_mc) {
        message("--- STARTING PARALLEL SIMULATION OF ", nchains, " CHAINS ---")
      } else {
        message("Performing parallel simulation of ", nchains, " chains...")
      }
    } else {
      if (.Platform$OS.type != "windows") devout <- '/dev/null' else devout <- 'nul:'
    }

    if (have_mc) {
      res <- parallel::mclapply(
        X = seq_len(nchains),
        FUN = bayesmr_mix_het_fit_parallel,
        mc.cores = threads,
        mc.set.seed = TRUE,
        data.c = data,
        control.c = control,
        prior.c = prior,
        lib = .bayesmrEnv$path.to.me
      )
    } else {
      if (is.null(cl)) {
        cl_local <- parallel::makePSOCKcluster(
          rep("localhost", threads),
          outfile = if (verbose) "" else devout
        )
        on.exit(parallel::stopCluster(cl_local), add = TRUE)
        cl_to_use <- cl_local

        ## Seed the cluster explicitly for reproducible independent streams
        parallel::clusterSetRNGStream(cl_to_use, iseed = internal_seed)
      } else {
        warning("for user-supplied clusters, the RNG seed must be set externally.")
        cl_to_use <- cl
      }

      res <- parallel::parLapply(
        cl = cl_to_use,
        X = seq_len(nchains),
        fun = bayesmr_mix_het_fit_parallel,
        data.c = data,
        control.c = control,
        prior.c = prior,
        lib = .bayesmrEnv$path.to.me
      )
    }
  
    if (verbose) {
      if (.Platform$OS.type != "windows" && !have_mc) {
        message("--- END OF PARALLEL SIMULATION OF ", nchains, " CHAINS ---\n")
      }
    }
  } else {
    ## Serial execution: ordinary RNG is enough
    if (!is.null(seed)) {
      set.seed(internal_seed)
    }

    res <- vector("list", nchains)
    for (ch in seq_len(nchains)) {
      res[[ch]] <- bayesmr_mix_het_fit_serial(
        ch = ch,
        data.c = data,
        control.c = control,
        prior.c = prior
      )
    }
  }

  # final post-processing of all chains
  if (nchains > 1L && post_all) {
    # gamma.chain <- res[[1]]@gamma.chain
    # beta.chain <- res[[1]]@beta.chain
    # xi.chain <- res[[1]]@xi.chain
    # alpha.chain <- res[[1]]@alpha.chain
    # niter <- length(gamma.chain)
    # for (ch in 2:nchains) {
    #   gamma.chain <- abind::abind(gamma.chain, res[[ch]]@gamma.chain, along = 1)
    #   beta.chain <- c(beta.chain, res[[ch]]@beta.chain)
    #   xi.chain <- abind::abind(xi.chain, res[[ch]]@xi.chain, along = 1)
    #   alpha.chain <- abind::abind(alpha.chain, res[[ch]]@alpha.chain, along = 1)
    # }

    # if (control[["verbose"]]) message("Final post-processing of all chains:")

    # if (control[["verbose"]]) {
    #   # message("done!")
    #   close(pb)
    # }

    # for (ch in 1:nchains) {
    #   res[[ch]]@gamma.chain <- gamma.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
    #   res[[ch]]@beta.chain <- beta.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
    #   res[[ch]]@xi.chain <- xi.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
    #   res[[ch]]@alpha.chain <- alpha.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
    # }
  }

  res <- new("bayesmr_mix_fit_list", results = res)

  return(res)
}
