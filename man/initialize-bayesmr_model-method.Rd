% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesmrclus_classes.R
\name{initialize,bayesmr_model-method}
\alias{initialize,bayesmr_model-method}
\alias{bayesmr_model-initialize}
\title{Create an instance of the \code{bayesmr_model} class using new/initialize.}
\usage{
\S4method{initialize}{bayesmr_model}(.Object, p = numeric(), G = numeric())
}
\arguments{
\item{.Object}{Prototype object from the class \code{\link{bayesmr_model}}.}

\item{p}{A length-one character vector representing the number of dimensions
of the latent space to use in the MDS analysis.}

\item{G}{A length-one numeric vector representing the number of clusters to
partition the subjects into.}
}
\description{
Create an instance of the \code{bayesmr_model} class using new/initialize.
}
\author{
Sergio Venturini \email{sergio.venturini@unicatt.it}
}
