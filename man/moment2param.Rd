% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moment2param.R
\name{moment2param}
\alias{moment2param}
\title{Moment to Parameters: A Batch Process}
\usage{
moment2param(distname, FUN = paste0("moment2", distname), ...)
}
\arguments{
\item{distname}{\link[base]{character} scalar, distribution name.
Currently supported are \code{'GH'} for Tukey \eqn{g}-&-\eqn{h} distribution,
\code{'sn'} for skew-normal distribution and \code{'st'} for skew-\eqn{t} distribution}

\item{FUN}{\link[base]{name} or \link[base]{character} scalar,
(name of) \link[base]{function} used to solve the distribution parameters from moments.
Default is \code{paste0('moment2', distname)}, e.g., \link{moment2GH} will be used for \code{distname = 'GH'}.
To use one of the educational functions, specify
\code{FUN = moment2GH_g_demo} or \code{FUN = 'moment2GH_g_demo'}.}

\item{...}{\link[base]{numeric} scalars,
some or all of \code{mean}, \code{sd}, \code{skewness} and \code{kurtosis}
(length will be recycled).}
}
\value{
Function \link{moment2param} returns a \link[base]{list} of \link[base]{numeric} \link[base]{vector}s.
}
\description{
Converts multiple sets of moments to multiple sets of distribution parameters.
}
\examples{
skw = c(.2, .5, .8)
krt = c(.5, 1, 1.5)
moment2param(distname = 'GH', skewness = skw, kurtosis = krt)
moment2param(distname = 'st', skewness = skw, kurtosis = krt)

}
