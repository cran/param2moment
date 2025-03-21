% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TukeyGH.R
\name{moment2GH}
\alias{moment2GH}
\alias{moment2GH_h_demo}
\alias{moment2GH_g_demo}
\title{Solve Tukey \eqn{g}-&-\eqn{h} Parameters from Moments}
\usage{
moment2GH(mean = 0, sd = 1, skewness, kurtosis)

moment2GH_h_demo(sd = 1, kurtosis)

moment2GH_g_demo(mean = 0, sd = 1, skewness)
}
\arguments{
\item{mean}{\link[base]{numeric} scalar, mean \eqn{\mu}, default value 0}

\item{sd}{\link[base]{numeric} scalar, standard deviation \eqn{\sigma}, default value 1}

\item{skewness}{\link[base]{numeric} scalar}

\item{kurtosis}{\link[base]{numeric} scalar}
}
\value{
Function \code{\link[=moment2GH]{moment2GH()}} returns a \link[base]{length}-4
\link[base]{numeric} \link[base]{vector} \eqn{(A, B, g, h)}.

Function \code{\link[=moment2GH_h_demo]{moment2GH_h_demo()}} returns a \link[base]{length}-2
\link[base]{numeric} \link[base]{vector} \eqn{(B, h)}.

Function \code{\link[=moment2GH_g_demo]{moment2GH_g_demo()}} returns a \link[base]{length}-3
\link[base]{numeric} \link[base]{vector} \eqn{(A, B, g)}.
}
\description{
Solve Tukey \eqn{g}-, \eqn{h}- and \eqn{g}-&-\eqn{h} distribution parameters
from mean, standard deviation, skewness and kurtosis.
}
\details{
Function \code{\link[=moment2GH]{moment2GH()}} solves the
location \eqn{A}, scale \eqn{B}, skewness \eqn{g}
and elongation \eqn{h} parameters of Tukey \eqn{g}-&-\eqn{h} distribution,
from user-specified mean \eqn{\mu} (default 0), standard deviation \eqn{\sigma} (default 1),
skewness and kurtosis.

An educational and demonstration function \code{\link[=moment2GH_h_demo]{moment2GH_h_demo()}} solves
\eqn{(B, h)} parameters of Tukey \eqn{h}-distribution,
from user-specified \eqn{\sigma} and kurtosis.
This is a non-skewed distribution, thus
the location parameter \eqn{A=\mu=0}, and the skewness parameter \eqn{g=0}.

An educational and demonstration function \code{\link[=moment2GH_g_demo]{moment2GH_g_demo()}} solves
\eqn{(A, B, g)} parameters of Tukey \eqn{g}-distribution,
from user-specified \eqn{\mu}, \eqn{\sigma} and skewness.
For this distribution, the elongation parameter \eqn{h=0}.
}
\examples{
moment2GH(skewness = .2, kurtosis = .3)

moment2GH_h_demo(kurtosis = .3)

moment2GH_g_demo(skewness = .2)

}
