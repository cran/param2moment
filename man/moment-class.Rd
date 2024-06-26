% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moment.R
\docType{class}
\name{moment-class}
\alias{moment-class}
\title{Raw, Central and Standardized Moments, and other Distribution Characteristics}
\description{
Up to 4th raw \eqn{\text{E}(Y^n)}, \href{https://en.wikipedia.org/wiki/Central_moment}{central} \eqn{\text{E}[(Y-\mu)^n]} and
\href{https://en.wikipedia.org/wiki/Standardized_moment}{standardized moments} \eqn{\text{E}[(Y-\mu)^n/\sigma^n]} of the random variable
\deqn{Y = (X - \text{location})/\text{scale}}

Also, the mean, standard deviation, skewness and excess kurtosis of the random variable \eqn{X}.
}
\details{
For \eqn{Y = (X - \text{location})/\text{scale}},
let \eqn{\mu = \text{E}(Y)}, then, according to
\href{https://en.wikipedia.org/wiki/Binomial_theorem}{Binomial theorem},
the 2nd to 4th central moments of \eqn{Y} are,
\deqn{\text{E}[(Y-\mu)^2] = \text{E}(Y^2) - 2\mu \text{E}(Y) + \mu^2 = \text{E}(Y^2) - \mu^2}
\deqn{\text{E}[(Y-\mu)^3] = \text{E}(Y^3) - 3\mu \text{E}(Y^2) + 3\mu^2 \text{E}(Y) - \mu^3 = \text{E}(Y^3) - 3\mu \text{E}(Y^2) + 2\mu^3}
\deqn{\text{E}[(Y-\mu)^4] = \text{E}(Y^4) - 4\mu \text{E}(Y^3) + 6\mu^2 \text{E}(Y^2) - 4\mu^3 \text{E}(Y) + \mu^4 = \text{E}(Y^4) - 4\mu \text{E}(Y^3) + 6\mu^2 \text{E}(Y^2) - 3\mu^4}

The distribution characteristics of \eqn{Y} are,
\deqn{\mu_Y = \mu}
\deqn{\sigma_Y = \sqrt{\text{E}[(Y-\mu)^2]}}
\deqn{\text{skewness}_Y = \text{E}[(Y-\mu)^3] / \sigma^3_Y}
\deqn{\text{kurtosis}_Y = \text{E}[(Y-\mu)^4] / \sigma^4_Y - 3}

The distribution characteristics of \eqn{X} are
\eqn{\mu_X = \text{location} + \text{scale}\cdot \mu_Y},
\eqn{\sigma_X = \text{scale}\cdot \sigma_Y},
\eqn{\text{skewness}_X = \text{skewness}_Y}, and
\eqn{\text{kurtosis}_X = \text{kurtosis}_Y}.
}
\section{Slots}{

\describe{
\item{\code{distname}}{\link[base]{character} scalar, name of distribution,
e.g., \code{'norm'} for normal, \code{'sn'} for skew-normal, \code{'st'} for skew-\eqn{t},
and \code{'GH'} for Tukey \eqn{g}-&-\eqn{h} distribution,
following the nomenclature of \link[stats]{dnorm}, \link[sn]{dsn}, \link[sn]{dst} and \code{QuantileGH::dGH}}

\item{\code{location,scale}}{\link[base]{numeric} scalars or \link[base]{vector}s,
location and scale parameters}

\item{\code{mu}}{\link[base]{numeric} scalar or \link[base]{vector},
1st \emph{raw} moment \eqn{\mu = \text{E}(Y)}.
Note that the 1st central moment \eqn{\text{E}(Y-\mu)} and
standardized moment \eqn{\text{E}(Y-\mu)/\sigma} are both 0.}

\item{\code{raw2,raw3,raw4}}{\link[base]{numeric} scalars or \link[base]{vector}s,
2nd or higher \emph{raw} moments \eqn{\text{E}(Y^n)}, \eqn{n\geq 2}}

\item{\code{central2,central3,central4}}{\link[base]{numeric} scalars or \link[base]{vector}s,
2nd or higher \href{https://en.wikipedia.org/wiki/Central_moment}{\emph{central} moments}, \eqn{\sigma^2 = \text{E}[(Y-\mu)^2]} and
\eqn{\text{E}[(Y-\mu)^n]}, \eqn{n\geq 3}}

\item{\code{standardized3,standardized4}}{\link[base]{numeric} scalars or \link[base]{vector}s,
3rd or higher \href{https://en.wikipedia.org/wiki/Standardized_moment}{\emph{standardized} moments},
\href{https://en.wikipedia.org/wiki/Skewness}{skewness} \eqn{\text{E}[(Y-\mu)^3]/\sigma^3} and
\href{https://en.wikipedia.org/wiki/Kurtosis}{kurtosis} \eqn{\text{E}[(Y-\mu)^4]/\sigma^4}.
Note that the 2nd standardized moment is 1}
}}

\note{
Potential name clash with function \code{e1071::moment}.
}
