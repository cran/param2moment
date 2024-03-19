

#' @title Moments of Skew-Normal Distribution
#' 
#' @description
#' Moments of \href{https://en.wikipedia.org/wiki/Skew_normal_distribution}{skew-normal distribution}, parameter nomenclature follows
#' \link[sn]{dsn} function.
#' 
#' @param xi,omega,alpha \link[base]{numeric} scalars or \link[base]{vector}s, 
#' location, scale and slant parameters of skew-normal distribution
#' 
#' @returns
#' Function [moment_sn] returns a \linkS4class{moment} object.
#' 
#' @examples
#' xi = 2; omega = 1.3; alpha = 3
#' moment_sn(xi, omega, alpha)
#' curve(sn::dsn(x, xi = 2, omega = 1.3, alpha = 3), from = 0, to = 6)
#' 
#' @importFrom methods new
#' @export
moment_sn <- function(xi = 0, omega = 1, alpha = 0) {
  do.call(what = new, args = c(list(Class = 'moment', location = xi, scale = omega), moment_sn_(alpha = alpha)))
}

moment_sn_ <- function(alpha = 0) {
  delta <- alpha / sqrt(1 + alpha^2)
  b <- sqrt(2/pi)
  mu <- b * delta
  moment_int(
    distname = 'sn', 
    mu = mu,
    raw2 = 1,
    raw3 = - pi/2 *mu^3 + 3*mu,
    raw4 = 3)
}



#' @title Solve Skew-Normal Parameters from Moments
#' 
#' @description
#' Solve skew-normal parameters from mean, standard deviation and skewness.
#' 
#' @param mean \link[base]{numeric} scalar, mean, default value 0
#' 
#' @param sd \link[base]{numeric} scalar, standard deviation, default value 1
#' 
#' @param skewness \link[base]{numeric} scalar
#' 
#' @returns
#' Function [param_sn] returns a \link[base]{numeric} \link[base]{vector} of 
#' \link[base]{length}-3, representing the 
#' location \eqn{\xi}, scale \eqn{\omega} and slant \eqn{\alpha} parameters.
#' 
#' @examples
#' param_sn(skewness = .3)
#' 
#' @importFrom stats optim
#' @export
param_sn <- function(mean = 0, sd = 1, skewness) {
  optim(par = c(xi = 0, omega = 1, alpha = 0), fn = function(x) {
    mm <- moment_sn_(alpha = x[3L])
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm)) - c(mean, sd, skewness))
  })$par
}


