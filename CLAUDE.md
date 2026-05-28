# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package overview

**BayesMRclus** is an R package (v0.0.7-7) for Bayesian inference in two-sample summary-data Mendelian Randomization (MR) with clustering of genetic variants (SNPs used as instrumental variables). It uses Rcpp/RcppArmadillo for the MCMC core (C++14, OpenMP) and S4 classes throughout.

## Build and development commands

All run from an R session in the project root:

```r
# Install dependencies and build
devtools::install_deps()
devtools::build()
devtools::install()

# Check (CRAN-style)
devtools::check()

# Generate documentation from roxygen2 comments
devtools::document()

# Compile C++ sources only (faster than full install during C++ development)
devtools::compile_dll()

# Load the package without installing (preferred during development)
devtools::load_all()
```

There is no formal testthat test suite. The `test/` directory contains manual example scripts (e.g., `test/examples.R`, `test/examples_mix_het/`).

## Architecture

### Model variants

There are four model variants, each with a corresponding top-level function:

| Function | Description |
|---|---|
| `bayesmr()` | No clustering, fixed heterogeneity |
| `bayesmr_het()` | No clustering, with heterogeneity (psi/tau parameters) |
| `bayesmr_mix()` | Mixture/clustering model, fixed heterogeneity |
| `bayesmr_mix_het()` | Mixture/clustering + heterogeneity (full model) |

All four share the same call pattern: `fn(data, control, prior, cl, post_all)`.

### Data flow

1. Input data is a data frame with columns `beta_exposure`, `beta_outcome`, `se_exposure`, `se_outcome`, `SNP`. The `bayesmr_data` S4 class (initialized via `new("bayesmr_data", data, n)`) validates, reorients alleles, and stores this.
2. `bayesmr_control()` and `bayesmr_prior()` produce named lists of MCMC and prior settings.
3. `bayesmr_init()` generates starting values.
4. The inner `bayesmr_fit()` / `bayesmr_het_fit()` etc. functions call the C++ MCMC routines (in `src/`) via `.Call()` and return a single-chain S4 fit object.
5. The outer `bayesmr()` / `bayesmr_mix_het()` etc. functions run multiple chains (optionally in parallel via `parallel::mclapply` or `parallel::parLapply`) and bundle results into a `*_fit_list` S4 object.

### S4 class hierarchy

- `bayesmr_data` — holds SNP summary statistics (slot `@data`: data frame, `@n`: number of SNPs)
- `bayesmr_fit` / `bayesmr_het_fit` / `bayesmr_mix_fit` / `bayesmr_mix_het_fit` — single chain results (slots: `@gamma.chain`, `@beta.chain`, `@xi.chain`, `@alpha.chain`, etc.)
- `bayesmr_fit_list` / `bayesmr_het_fit_list` / `bayesmr_mix_fit_list` / `bayesmr_mix_het_fit_list` — multi-chain wrappers (slot `@results`: list of single-chain fit objects)

Standard S4 generics (`show`, `summary`, `plot`, `[`) are implemented for all classes.

### R source files (`R/`)

| File | Responsibility |
|---|---|
| `bayesmrclus_classes.R` | All S4 class definitions and methods |
| `bayesmrclus.R` | `bayesmr()` — outer loop, parallelism |
| `bayesmrclus_het.R` | `bayesmr_het()` — outer loop |
| `bayesmrclus_mix.R` | `bayesmr_mix()` — outer loop |
| `bayesmrclus_mix_het.R` | `bayesmr_mix_het()` — outer loop |
| `bayesmrclus_fit.R` | Inner single-chain fitting functions |
| `bayesmrclus_control.R` | `bayesmr_control()`, `check_control()` |
| `bayesmrclus_prior.R` | `bayesmr_prior()`, `check_prior()` |
| `bayesmrclus_init.R` | `bayesmr_init()` — starting values |
| `bayesmrclus_lists.R` | `*_fit_list` class definitions |
| `bnpmr_postprocess.R` | Post-processing utilities (`read_mcmc_output`, `salso_partition`, etc.) |
| `cluster_agreement.R` | Cluster comparison metrics (ARI, VI, FMI, PSM, alluvial plots) |
| `mrclust_results.R` | Interface with the `mrclust` package for comparison |
| `estim.R`, `jointpost.R` | Likelihood and posterior computations in R |
| `distrib.R` | Half-t distribution functions (`dhalft`, `rhalft`, `EXhalft`) |

### C++ source files (`src/`)

The MCMC samplers are implemented in C++ using Rcpp (no Rcpp attributes — functions are exported via `RcppExport SEXP` and registered in `src/init.cpp`):

| File | Implements |
|---|---|
| `bayesmr.cpp` | No-clustering sampler |
| `mcmc.cpp` | Clustering sampler (fixed heterogeneity) |
| `mcmc_het.cpp` | No-clustering sampler with heterogeneity |
| `mcmc_mix.cpp` | Mixture sampler |
| `mcmc_mix_het.cpp` | Mixture + heterogeneity sampler (full model) |
| `mcmc_mix_het_2.cpp` | Alternative variant |
| `bayesmr_post.cpp` | Posterior evaluation utilities |
| `DPmix.cpp` | Dirichlet Process mixture utilities |
| `distrib.cpp` | C++ distribution helpers |
| `matrix_utils.cpp` | Matrix operations |
| `utils.cpp` | General utilities |
| `init.cpp` | Rcpp function registration table |

The header `inst/include/bayesmr.h` declares all exported C++ functions and is included by all `.cpp` files.

### Key parameters

`bayesmr_control()` defaults: `nsim=5000, burnin=10000, thin=1, nchains=1, threads=1, parallel="no", beta.prop=0.5, psi.prop=0.1, tau.prop=0.1`

`bayesmr_prior()` defaults: half-t priors on psi/tau (`alpha=0.1, nu=3`), normal priors on gamma/beta (`mean=0, var=1`), Beta(2,2) on alpha.

### Datasets

Real GWAS summary-data pairs are bundled in `data/` as `.rda` files (e.g., `ldl_cad`, `bmi_t2d`, `hdl_chd`, `whradjbmi_insulin`). Each is a data frame with SNP-level summary statistics. Example scripts in `test/data/` show how each dataset was prepared.
