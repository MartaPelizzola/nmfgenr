#' Unit deviance from the Tweedie distribution
#'
#' @param y A vector of observations.
#' @param mu A vector of estimates.
#' @param pwr Power parameter of the Tweedie distribution.
#'
#' @returns The unit deviance value from the Tweedie distribution between observations and estimates.
#' @export
#'
#' @examples
#' y <- rnorm(20,5)
#' mu <- rep(5,20)
#' unitdev(y,mu,0)
unitdev <- function(y, mu, pwr){
  if (length(mu)!=length(y)){
    stop("y and mu need to have the same length.")
  }
  r <- mu
  p <- which(y > 0)
  if (pwr == 2){
    r[p] <- y[p]/mu[p] - log(y[p]) + log(mu[p]) -1
  }
  else if(pwr == 1){
    r[p] <- (y * ( log(y) - log(mu)) - y + mu)[p]
  }
  else{
    r[p] <- ((y^(2-pwr) + (1 - pwr)*mu^(2-pwr) - (2-pwr)*y*mu^(1 -pwr))/((2-pwr)*(1-pwr)))[p]
  }
  return(sum(r))
}

#' Negative Binomial divergence.
#'
#' @param y A vector of observations.
#' @param mu A vector of estimates.
#' @param alpha_n A vector of dispersion parameter (one for each observation).
#'
#' @returns The Negative Binomial divergence between a vector of observations y and one of estimates mu.
#' @export
#'
#' @examples
#'
nbdev <- function(y, mu, alpha_n){
  if (length(mu)!=length(y) | length(alpha_n)!=length(y)){
    stop("y, mu and alpha_n need to have the same length.")
  }
  r <- mu
  p <- which(y > 0)
  r[p] <- (y * (log(y)- log(mu)) - (alpha_n + y) * (log(alpha_n + y) - log(alpha_n + mu) ) )[p]
  return(sum(r))
}


#' Generalized Kullback-Leibler divergence.
#'
#' @param y A vector of observations.
#' @param mu A vector of estimates.
#'
#' @returns The Generalized Kullback-Leibler divergence between a vector of observations y and one of estimates mu.
#' @export
#'
#' @examples
gkldev <- function(y, mu){
  r <- mu
  p <- which(y > 0)
  r[p] <- (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}
