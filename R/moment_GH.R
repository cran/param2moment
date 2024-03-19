



#' @title Moments of Tukey \eqn{g}-&-\eqn{h} Distribution
#' 
#' @description
#' Moments of Tukey \eqn{g}-&-\eqn{h} distribution.
#' 
#' @param A,B,g,h \link[base]{numeric} scalars or \link[base]{vector}s,
#' location, scale, skewness and elongation parameters of
#' Tukey \eqn{g}-&-\eqn{h} distribution 
#' 
#' @returns
#' Function [moment_GH] returns a \linkS4class{moment} object.
#' 
#' @references 
#' Raw moments of Tukey \eqn{g}-&-\eqn{h} distribution: \doi{10.1002/9781118150702.ch11}
#' 
#' @examples
#' A = 3; B = 1.5; g = .7; h = .01
#' moment_GH(A = A, B = B, g = 0, h = h)
#' moment_GH(A = A, B = B, g = g, h = 0)
#' moment_GH(A = A, B = B, g = g, h = h)
#' 
#' @importFrom methods new
#' @export
moment_GH <- function(A = 0, B = 1, g = 0, h = 0) {
  do.call(what = new, args = c(list(Class = 'moment', location = A, scale = B), moment_GH_(g = g, h = h)))
}

moment_GH_ <- function(g = 0, h = 0) {
  tmp <- data.frame(g, h) # recycling
  g <- tmp[[1L]]
  h <- tmp[[2L]]
  
  g0 <- (g == 0)
  
  # 1st-4th raw moment E(Y^n), when `g = 0`
  mu <- r3 <- rep(0, times = length(g))
  r2 <- 1 / (1-2*h) ^ (3/2) # (45a), page 502
  r4 <- 3 / (1-4*h) ^ (5/2) # (45b), page 502
  
  if (any(!g0)) {
    # {r}aw moment E(Y^n), when `g != 0`
    r_g <- function(n) { # equation (47), page 503
      tmp <- lapply(0:n, FUN = function(i) {
        (-1)^i * choose(n,i) * exp((n-i)^2 * g^2 / 2 / (1-n*h))
      })
      suppressWarnings(Reduce(f = `+`, tmp) / g^n / sqrt(1-n*h)) # warnings for `h > 1/n`
    }
    
    mu[!g0] <- r_g(1L)[!g0]
    r2[!g0] <- r_g(2L)[!g0]
    r3[!g0] <- r_g(3L)[!g0]
    r4[!g0] <- r_g(4L)[!g0]
  }
  
  moment_int(distname = 'GH', mu = mu, raw2 = r2, raw3 = r3, raw4 = r4)
  
}



#' @title Solve Tukey \eqn{g}-&-\eqn{h} Parameters from Moments
#' 
#' @description
#' Solve Tukey \eqn{g}-, \eqn{h}- and \eqn{g}-&-\eqn{h} distribution parameters 
#' from mean, standard deviation, skewness and kurtosis.
#' 
#' @param mean \link[base]{numeric} scalar, mean, default value 0
#' 
#' @param sd \link[base]{numeric} scalar, standard deviation, default value 1
#' 
#' @param skewness \link[base]{numeric} scalar
#' 
#' @param kurtosis \link[base]{numeric} scalar
#' 
#' @details
#' Function [param_GH] solves the 
#' location \eqn{A}, scale \eqn{B}, skewness \eqn{g} 
#' and elongation \eqn{h} parameters of Tukey \eqn{g}-&-\eqn{h} distribution,
#' from user-specified mean \eqn{\mu}, standard deviation \eqn{\sigma}, 
#' skewness and kurtosis.  
#' 
#' @returns
#' Function [param_GH] returns a \link[base]{length}-4 
#' \link[base]{numeric} \link[base]{vector} \eqn{(A, B, g, h)}.
#' 
#' @examples
#' param_GH(skewness = .2, kurtosis = .3)
#' 
#' @name param_GH
#' @importFrom stats optim
#' @export
param_GH <- function(mean = 0, sd = 1, skewness, kurtosis) {
  optim(par = c(A = 0, B = 1, g = 0, h = 0), fn = function(x) {
    mm <- moment_GH_(g = x[3L], h = x[4L])
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm), kurtosis_moment_(mm)) - c(mean, sd, skewness, kurtosis))
  })$par
}


#' @rdname param_GH
#' 
#' @details
#' A education and demonstration function [param_GH_h] solves the 
#' scale \eqn{B} and elongation \eqn{h} parameters of Tukey \eqn{h}-distribution,
#' from user-specified standard deviation \eqn{\sigma} and kurtosis.
#' This is a non-skewed distribution, thus 
#' the location parameter \eqn{A} is the mean \eqn{\mu}, and the skewness parameter \eqn{g=0}.
#' 
#' @returns
#' Function [param_GH_h] returns a \link[base]{length}-2 
#' \link[base]{numeric} \link[base]{vector} \eqn{(B, h)}.
#' 
#' @examples
#' param_GH_h(kurtosis = .3)
#' 
#' @importFrom stats optim
#' @export
param_GH_h <- function(sd = 1, kurtosis) {
  optim(par = c(B = 1, h = 0), fn = function(x) {
    mm <- moment_GH_(g = 0, h = x[2L])
    crossprod(c(sd_moment_(mm, scale = x[1L]), kurtosis_moment_(mm)) - c(sd, kurtosis))
  })$par
}

#' @rdname param_GH
#' 
#' @details
#' A education and demonstration function [param_GH_g] solves the 
#' location \eqn{A}, scale \eqn{B} and skewness \eqn{g} parameters of Tukey \eqn{g}-distribution,
#' from user-specified mean \eqn{\mu}, standard deviation \eqn{\sigma} and skewness.
#' For this distribution, the elongation parameter \eqn{h=0}.
#' 
#' @returns
#' Function [param_GH_g] returns a \link[base]{length}-3 
#' \link[base]{numeric} \link[base]{vector} \eqn{(A, B, g)}.
#' 
#' @examples
#' param_GH_g(skewness = .2)
#' 
#' @importFrom stats optim
#' @export
param_GH_g <- function(mean = 0, sd = 1, skewness) {
  optim(par = c(A = 0, B = 1, g = 0), fn = function(x) {
    mm <- moment_GH_(g = x[3L], h = 0)
    crossprod(c(mean_moment_(mm, location = x[1L], scale = x[2L]), sd_moment_(mm, scale = x[2L]), skewness_moment_(mm)) - c(mean, sd, skewness))
  })$par
}



