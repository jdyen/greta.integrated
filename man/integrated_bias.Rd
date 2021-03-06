% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/integrated_bias-class.R
\name{integrated_bias}
\alias{integrated_bias}
\alias{no_bias}
\alias{detection}
\alias{is.integrated_bias}
\alias{print.integrated_bias}
\alias{summary.integrated_bias}
\title{create integrated bias objects}
\usage{
no_bias()

detection(detection = beta(1, 1))

is.integrated_bias(object)

\method{print}{integrated_bias}(x, ...)

\method{summary}{integrated_bias}(object, ...)
}
\arguments{
\item{object}{an \code{integrated_bias} object}

\item{x}{}

\item{...}{additional arguments to \link[base]{print}, \link[base]{summary}, and \link[graphics]{plot} methods (currently ignored)}

\item{p_detect}{parameter of detection model (see details for information on setting prior distributions)}

\item{x}{an \code{integrated_bias} object}
}
\value{
An object of class \code{integrated_bias}, which can be used to create
   \link[greta.integrated]{integrated_data} and \link[greta.integrated]{integrated_model} objects
}
\description{
An \code{integrated_bias} object contains a bias model
  that can be shared among multiple components in an integrated population
  model
}
\details{
Prior distributions can be specified as single-dimensional
  greta distribution, e.g., \code{normal(0, 1)}. Link functions and transformations
  can be specified directly in-line, e.g., \code{ilogit(normal(0, 1))} specifies
  normal priors with a mean of zero and a standard deviation of one, transformed
  with an inverse-logit link.
}
\examples{
\dontrun{

library(integrated)

# a really basic age-structured model with five age classes
process <- leslie(5, density = "none")

# setting custom priors
process <- detection(params = list(abundance = ilogit(normal(0, 1))))
}
}
