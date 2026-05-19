#' Estimate a BayesMR model.
#'
#' `bayesmr()` is the main fitting function of the package. It estimates a
#' BayesMR model from an object of class [`bayesmr_data-class`] using MCMC.
#' Multiple chains can be run serially or in parallel according to the settings
#' supplied through [`bayesmr_control()`].
#'
#' @param data An object of class [`bayesmr_data-class`] containing the summary
#'   data to be analyzed.
#' @param control A named list of control parameters for the MCMC sampler, such
#'   as the number of iterations, burn-in, thinning, number of chains, seed, and
#'   parallelization options. Missing entries are filled with the defaults from
#'   [`bayesmr_control()`].
#' @param prior A named list of prior hyperparameters, or `NULL`. If `NULL`, the
#'   default prior specification returned by [`bayesmr_prior()`] is used. Missing
#'   entries are filled with the defaults from [`bayesmr_prior()`].
#' @param cl An optional cluster object created by [parallel::makeCluster()] or
#'   a compatible `snow` cluster. This is used only when
#'   `control$parallel = "snow"` and more than one chain is run. If `NULL`, a
#'   temporary local PSOCK cluster is created and stopped when the function exits.
#'   For user-supplied clusters, the random-number seed must be set externally.
#' @param post_all A length-one logical vector. Currently reserved for optional
#'   post-processing of multiple simulated chains.
#'
#' @return An object of class [`bayesmr_fit_list-class`] containing one fitted
#'   [`bayesmr_fit-class`] object for each simulated chain.
#'
#' @details
#' The function validates the supplied `control` and `prior` lists, initializes
#' each chain with [`bayesmr_init()`], and fits each chain with `bayesmr_fit()`.
#' If a seed is supplied through `control$seed`, it is used internally while
#' preserving the caller's random-number generator state.
#'
#' Parallel execution is used only when more than one chain is requested and
#' `control$parallel` is set to either `"multicore"` or `"snow"`. Multicore
#' execution is unavailable on Windows.
#'
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#'
#' @seealso
#' [`bayesmr_data-class`] for the input data format;
#' [`bayesmr_fit_list-class`] for the returned object;
#' [`bayesmr_control()`] and [`bayesmr_prior()`] for control and prior settings.
#'
#' @references
#' Consonni, G., Venturini, S., Castelletti, F. (2026).
#' "Bayesian Hierarchical Modeling for Two-Sample Summary-Data Mendelian
#' Randomization under Heterogeneity". Working paper.
#'
#' @examples
#' \dontrun{
#' data(ldl_cad)
#' dat <- new("bayesmr_data", ldl_cad)
#'
#' control <- bayesmr_control(
#'   nsim = 1000,
#'   burnin = 1000,
#'   nchains = 2,
#'   verbose = TRUE
#' )
#' prior <- bayesmr_prior()
#'
#' res <- bayesmr(dat, control = control, prior = prior)
#' summary(res)
#' plot(res, what = "trace")
#' }
#'
#' @export
bayesmr <- function(data, control = bayesmr_control(), prior = NULL, cl = NULL, post_all = FALSE) {
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
  if (is.null(control[["K_start"]]) || control[["K_start"]] > 1) {
    control[["K_start"]] <- 1
    message("initial number of clusters set to one.")
  }
  K_start <- control[["K_start"]]
  store.burnin <- control[["store.burnin"]]
  verbose <- control[["verbose"]]

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

  bayesmr_fit_parallel <- function(c, data.c, control.c, prior.c, lib) {
    suppressMessages(require(BayesMRclus, lib.loc = lib))
    control.c[["verbose"]] <- FALSE

    start.c <- bayesmr_init(
      data = data.c,
      random_start = control.c[["random_start"]],
      K_start = control.c[["K_start"]]
    )
    
    bayesmr_fit(
      data = data.c,
      control = control.c,
      prior = prior.c,
      start = start.c
    )
  }

  bayesmr_fit_serial <- function(ch, data.c, control.c, prior.c) {
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
    
    fit <- bayesmr_fit(
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
        FUN = bayesmr_fit_parallel,
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
        fun = bayesmr_fit_parallel,
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
      res[[ch]] <- bayesmr_fit_serial(
        ch = ch,
        data.c = data,
        control.c = control,
        prior.c = prior
      )
    }
  }

  res <- new("bayesmr_fit_list", results = res)

  return(res)
}
