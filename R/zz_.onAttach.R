# Package-wide global variables
.bayesmrEnv <- new.env()

#' @importFrom utils globalVariables
#' @importFrom tools file_path_as_absolute
.onLoad <- function(lib, pkg){
  # needed to avoid annoying notes in R CMD CHECK
  # (see https://github.com/tidyverse/magrittr/issues/29)
  if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(".", ""))
  }
  .bayesmrEnv$path.to.me <- tools::file_path_as_absolute(lib)
  .bayesmrEnv$nlog.double.eps <- -log(.Machine[["double.eps"]])
  .bayesmrEnv$current_p <- 2
  .bayesmrEnv$current_G <- 3
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(sprintf("Package %s (%s) loaded.\nTo cite, type citation(\"%s\")",
    pkgname, utils::packageDescription(pkgname)$Version, pkgname))
  # packageStartupMessage("Caution: ", pkgname, " is under active developement! Breaking changes may occur in the future.")
  packageStartupMessage("Please report improvements and bugs to: https://github.com/sergioventurini/BayesMRclus/issues")
}

.onUnload <- function(lib) {
  library.dynam.unload("BayesMRclus", lib)
}
