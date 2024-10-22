% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayesmrclus_classes.R
\name{plot,bayesmr_data-method}
\alias{plot,bayesmr_data-method}
\alias{bayesmr_data-plot}
\title{Provide a graphical summary of a \code{bayesmr_data} class instance.}
\usage{
\S4method{plot}{bayesmr_data}(x, colors = c("white", "black"), font = NA, cex.font = NA, ...)
}
\arguments{
\item{x}{An object of class \code{\link{bayesmr_data}}.}

\item{colors}{A character vector providing the colors to use in the plot.}

\item{font}{A length-one numeric vector for the font to use for text.
Can be a vector. \code{NA} values (the default) mean use \code{par("font")}.}

\item{cex.font}{A length-one numeric vector for the character expansion
factor. \code{NULL} and \code{NA} are equivalent to \code{1.0}. This is an
absolute measure, not scaled by \code{par("cex")} or by setting
'   \code{par("mfrow")} or \code{par("mfcol")}. Can be a vector.}

\item{...}{Further arguments to pass on (currently ignored).}
}
\description{
Provide a graphical summary of a \code{bayesmr_data} class instance.
}
\examples{
data(simdiss)
library(bayesplot)
cols <- color_scheme_set("brightblue")
plot(simdiss, colors = unlist(cols)[c(1, 6)], font = 1, cex.font = 0.75)
}
\author{
Sergio Venturini \email{sergio.venturini@unicatt.it}
}
