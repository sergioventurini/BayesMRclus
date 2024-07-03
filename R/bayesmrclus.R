#' Estimation of a BayesMR model.
#'
#' \code{bayesmr()}, the main function of the package, estimates a BayesMR model
#'   for a given set of \emph{S} dissimilarity matrices.
#'
#' @param data An object of class \code{bayesmr_data} containing the data
#'   to analyze.
#' @param p A length-one numeric vector indicating the number of dimensions of the
#'   data space.
#' @param G A length-one numeric vector indicating the number of cluster to
#'   partition the \emph{S} subjects.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution. See
#'   \code{\link{bayesmr_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{bayesmr_prior}()} for more details.
#' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use if
#'   \code{parallel = "snow"}. If not supplied, a cluster on the local machine
#'   is created for the duration of the \code{bayesmr()} call.
#' @param post_all A length-one logical vector, which if TRUE applies a further
#'   post-processing to the simulated chains (in case these are more than one).
#' @return A \code{bayesmr_fit_list} object.
#' @author Sergio Venturini \email{sergio.venturini@unicatt.it}
#' @seealso \code{\link{bmds}} for Bayesian (metric) multidimensional scaling.
#' @seealso \code{\link{bayesmr_data}} for a description of the data format.
#' @seealso \code{\link{bayesmr_fit_list}} for a description of the elements
#'   included in the returned object.
#' @references
#'   [[TODO]]
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.bayesmr <- bayesmr(simdiss, p, G, control)
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
bayesmr <- function(data, p = 1, G = 1, control = bayesmr_control(), prior = NULL, cl = NULL, post_all = FALSE) {
  if (p < 1)
    stop("the number of data dimensions p must be at least one.")
  if (G < 1)
    stop("the number of clusters/groups G must be at least one.")
  
  .bayesmrEnv$current_p <- p
  .bayesmrEnv$current_G <- G

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
  random.start <- control[["random.start"]]
  partition <- control[["partition"]]
  store.burnin <- control[["store.burnin"]]
  verbose <- control[["verbose"]]

  have_mc <- have_snow <- FALSE
  if (parallel != "no" && threads > 1L) {
    if (parallel == "multicore") have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) {
      warning("number of cores forced to 1 (i.e. no parallel computing used).")
      threads <- 1L
    }
    loadNamespace("parallel") # get this out of the way before recording seed
  }

  n <- nrow(data)
  totiter <- burnin + nsim
  p <- as.integer(p)
  G <- as.integer(G)
  
  # save current random number generator kind
  old.rng <- RNGkind()[1L]
  RNGkind(kind = "L'Ecuyer-CMRG")
  
  # perform MCMC simulation
  if (nchains > 1L && (have_mc || have_snow)) {
    bayesmr_fit_parallel <- function(c, data.c, p.c, G.c, control.c, prior.c, lib) {
      suppressMessages(require(bayesmr, lib.loc = lib))
      control.c[["verbose"]] <- FALSE
      # message("Starting cluster node ", c, " on local machine")
      start.c <- bayesmr_init(data = data.c, p = p.c, G = G.c, random.start = control.c[["random.start"]],
        partition = control.c[["partition"]])
      bayesmr_fit(data = data.c, p = p.c, G = G.c, control = control.c, prior = prior.c, start = start.c)
    }
    # environment(bayesmr_fit_parallel) <- .GlobalEnv # this prevents passing objects other than those needed for
    #                                                 # evaluating the bayesmr_fit_parallel function

    if (is.null(prior)) {
      prior <- bayesmr_prior()
    } else {
      prior <- check_list_na(prior, bayesmr_prior())
    }
    if (!check_prior(prior)) {
      stop("the prior hyperparameter list is not correct; see the documentation for more details.")
    }

    if (verbose) {
      devout <- ""
      if (.Platform$OS.type != "windows" && !have_mc) {
        message("--- STARTING PARALLEL SIMULATION OF ", nchains, " CHAINS ---")
      } else {
        message("Performing parallel simulation of ", nchains, " chains...")
      }
    } else {
      if (.Platform$OS.type != "windows") {
        devout <- '/dev/null'
      } else {
        devout <- 'nul:'
      }
    }

    res <- if (have_mc) {
             if (!is.null(seed)) {
               set.seed(seed)
               parallel::mc.reset.stream()
             }
             parallel::mclapply(seq_len(nchains), bayesmr_fit_parallel, mc.cores = threads, mc.set.seed = TRUE,
               data.c = data, p.c = p, G.c = G, control.c = control, prior.c = prior, lib = .bayesmrEnv$path.to.me)
           } else if (have_snow) {
             if (is.null(cl)) {
               cl <- parallel::makePSOCKcluster(rep("localhost", threads), outfile = devout) # outfile doesn't work on 
                                                                                             # Windows
               parallel::clusterSetRNGStream(cl, seed)
               res <- parallel::parLapply(cl, seq_len(nchains), bayesmr_fit_parallel, data.c = data, p.c = p, G.c = G, 
                 control.c = control, prior.c = prior, lib = .bayesmrEnv$path.to.me)
               parallel::stopCluster(cl)
               res
             } else parallel::parLapply(cl, seq_len(nchains), bayesmr_fit_parallel, data.c = data, p.c = p, G.c = G,
               control.c = control, prior.c = prior, lib = .bayesmrEnv$path.to.me)
           }

    if (verbose) {
      if (.Platform$OS.type != "windows" && !have_mc){
        message("--- END OF PARALLEL SIMULATION OF ", nchains, " CHAINS ---\n")
      } else {
        # message("done!")
      }
    }
  } else {
    if (!is.null(seed)) {
     set.seed(seed)
    }
    res <- list()
    for (ch in 1:nchains) {
      if (verbose && nchains > 1L) message("--- STARTING SIMULATION OF CHAIN ", ch, " OF ", nchains, " ---")

      if (verbose) message("Initialization of the algorithm...")
  
      bayesmr.start <- bayesmr_init(data, p, G, random.start, partition = partition)
      if (is.null(prior)) {
        prior <- bayesmr_prior()
      } else {
        prior <- check_list_na(prior, bayesmr_prior())
      }
      if (!check_prior(prior))
        stop("the prior hyperparameter list is not correct; see the documentation for more details.")
    
      if (verbose) {
        # message("done!")
      }

      res[[ch]] <- bayesmr_fit(data = data, p = p, G = G, control = control, prior = prior, start = bayesmr.start)

      if (verbose && nchains > 1L) message("--- END OF CHAIN ", ch, " OF ", nchains, " ---\n")
    }
  }

  # final post-processing of all chains:
  if (nchains > 1 && post_all) {
    gamma.chain <- res[[1]]@gamma.chain
    beta.chain <- res[[1]]@beta.chain
    # z.chain <- res[[1]]@z.chain
    # z.chain.p <- res[[1]]@z.chain.p
    # alpha.chain <- res[[1]]@alpha.chain
    # eta.chain <- res[[1]]@eta.chain
    # sigma2.chain <- res[[1]]@sigma2.chain
    # lambda.chain <- res[[1]]@lambda.chain
    # prob.chain <- res[[1]]@prob.chain
    # x.ind.chain <- res[[1]]@x.ind.chain
    niter <- dim(z.chain.p)[1]
    for (ch in 2:nchains) {
      gamma.chain <- abind::abind(gamma.chain, res[[ch]]@gamma.chain, along = 1)
      beta.chain <- abind::abind(beta.chain, res[[ch]]@beta.chain, along = 1)
      # z.chain <- abind::abind(z.chain, res[[ch]]@z.chain, along = 1)
      # z.chain.p <- abind::abind(z.chain.p, res[[ch]]@z.chain.p, along = 1)
      # alpha.chain <- abind::abind(alpha.chain, res[[ch]]@alpha.chain, along = 1)
      # eta.chain <- abind::abind(eta.chain, res[[ch]]@eta.chain, along = 1)
      # sigma2.chain <- abind::abind(sigma2.chain, res[[ch]]@sigma2.chain, along = 1)
      # lambda.chain <- abind::abind(lambda.chain, res[[ch]]@lambda.chain, along = 1)
      # prob.chain <- abind::abind(prob.chain, res[[ch]]@prob.chain, along = 1)
      # x.ind.chain <- abind::abind(x.ind.chain, res[[ch]]@x.ind.chain, along = 1)
    }

    if (control[["verbose"]]) message("Final post-processing of all chains:")

    # ## Procrustes transformation of Z_g
    # if (control[["verbose"]]) message("   - applying Procrustes transformation...")
    # if (control[["verbose"]])
    #   pb <- bayesmr_pb(min = 0, max = (niter*G*nchains - 1), width = 49)
    # no <- 0
    # for (it in 1:(niter*nchains)) {
    #   for (g in 1:G) {
    #     if (control[["verbose"]]) bayesmr_setpb(pb, no)
    #     if (p == 1) {
    #       z.chain.p[it, , , g] <- as.numeric(procrustes(as.matrix(z.chain[it, , , g]),
    #         as.matrix(z.chain.p[(niter*nchains), , , g]), translation = TRUE, dilation = FALSE)$X.new)
    #     } else {
    #       z.chain.p[it, , , g] <- procrustes(z.chain[it, , , g], z.chain.p[(niter*nchains), , , g],
    #         translation = TRUE, dilation = FALSE)$X.new
    #     }
    #     no <- no + 1
    #   }
    # }
    if (control[["verbose"]]) {
      # message("done!")
      close(pb)
    }

    # if (G > 1) {
    #   if ((niter*nchains) > 10) {
    #     if (control[["verbose"]]) message("   - relabeling the parameter chain...")
    #     init <- ifelse((niter*nchains) <= 100, 5, 100)
        
    #     theta <- .Call('bayesmr_pack_par', PACKAGE = 'bayesmr',
    #       radz = as.double(z.chain.p),
    #       radalpha = as.double(alpha.chain),
    #       radlambda = as.double(lambda.chain),
    #       rn = as.integer(n),
    #       rp = as.integer(p),
    #       rM = as.integer((niter*nchains)),
    #       rG = as.integer(G)
    #     )

    #     theta.relab <- .Call('bayesmr_relabel', PACKAGE = 'bayesmr',
    #       radtheta = as.double(theta),
    #       radz = as.double(z.chain.p),
    #       radalpha = as.double(alpha.chain),
    #       radeta = as.double(eta.chain),
    #       radsigma2 = as.double(sigma2.chain),
    #       radlambda = as.double(lambda.chain),
    #       radprob = as.double(prob.chain),
    #       raix_ind = as.integer(x.ind.chain),
    #       rinit = as.integer(init),
    #       rn = as.integer(n),
    #       rp = as.integer(p),
    #       rS = as.integer(S),
    #       rM = as.integer((niter*nchains)),
    #       rR = as.integer(m + 1),
    #       rG = as.integer(G),
    #       rverbose = as.integer(control[["verbose"]])
    #     )

    #     theta <- array(theta.relab[[1]], c((niter*nchains), (m + 1), G))  # this is not needed elsewhere
    #     z.chain.p <- array(theta.relab[[2]], c((niter*nchains), n, p, G))
    #     alpha.chain <- array(theta.relab[[3]], c((niter*nchains), G))
    #     eta.chain <- array(theta.relab[[4]], c((niter*nchains), G))
    #     sigma2.chain <- array(theta.relab[[5]], c((niter*nchains), G))
    #     lambda.chain <- array(theta.relab[[6]], c((niter*nchains), G))
    #     prob.chain <- array(theta.relab[[7]], c((niter*nchains), S, G))
    #     x.ind.chain <- array(theta.relab[[8]], c((niter*nchains), S, G))
    #     x.chain <- t(apply(x.ind.chain, 1, function(x) as.integer(x %*% 1:G)))

    #     # if (control[["verbose"]]) # message("done!")
    #   } else {
    #     warning("the number of iterations is too small for relabeling; relabeling skipped.", call. = FALSE,
    #       immediate. = TRUE)
    #   }
    # }

    for (ch in 1:nchains) {
      res[[ch]]@gamma.chain <- gamma.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
      res[[ch]]@beta.chain <- beta.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
      # res[[ch]]@z.chain.p <- z.chain.p[(niter*(ch - 1) + 1):(niter*ch), , , , drop = FALSE]
      # res[[ch]]@alpha.chain <- alpha.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
      # res[[ch]]@eta.chain <- eta.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
      # res[[ch]]@sigma2.chain <- sigma2.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
      # res[[ch]]@lambda.chain <- lambda.chain[(niter*(ch - 1) + 1):(niter*ch), , drop = FALSE]
      # res[[ch]]@prob.chain <- prob.chain[(niter*(ch - 1) + 1):(niter*ch), , , drop = FALSE]
      # res[[ch]]@x.ind.chain <- x.ind.chain[(niter*(ch - 1) + 1):(niter*ch), , , drop = FALSE]
    }
  }

  # restore previous random number generator kind
  RNGkind(kind = old.rng)

  res <- new("bayesmr_fit_list", results = res)

  return(res)
}
