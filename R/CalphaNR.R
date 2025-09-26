#' @importFrom stats var
#' @importFrom stats runif
#'
#' @title Estimation of the dispersion parameter for Negative Binomial convex non-negative matrix factorization.
#'
#' @description Likelihood estimation of the dispersion parameter in Negative Binomial using Newton-Raphson.
#' The dispersion parameter can be either a vector of length one or vector of length no. of observation if obs-specific dispersion is used.
#'
#'
#' @param data Numeric matrix of counts data.
#' @param k Rank to be used for the convex non-negative matrix factorization
#' @param specific Logical. If TRUE obs-spcific dispersion is used in the Negative Binomial model.
#'
#'
#' @return Dispersion parameter for convex Negative Binomial NMF. Either vector of length one or vector of length no. of observations if obs-specific dispersion is used.
#'
#' @examples
#' alpha <- convex_alphaNR()
#'
#' @export
#'
#'
convex_alphaNR <- function(data, k=NULL, specific = FALSE){
  if (k!=round(k)){
    stop("The rank must be an integer.")
  }
  if(is.null(data)){
    stop("The data set is missing.")
  }
  if(is.null(k)){
    stop("A value for the rank is missing.")
  }
  if(length(k)!=1){
    stop("'k' has length larger than 1.")
  }
  res_p <- CNMFTweedie(t(data), K = k, pwr = 1)
  h_p <- t(res_p$E) %*% data
  w_p <- t(res_p$D)

  # differentiated once
  neglikdiff1 = function(alpha, data, estimate){
    sum(digamma(data + alpha) - digamma(alpha) - data/(alpha+estimate) - alpha/(alpha+estimate) + log(alpha/(alpha+estimate)) + 1)
  }

  # differentiated twice
  neglikdiff2 = function(alpha, data, estimate){
    sum(trigamma(data + alpha) - trigamma(alpha) + data/(alpha+estimate)^2 + 1/alpha - 2/(alpha+estimate) + alpha/(alpha+estimate)^2)
  }

  NR_alpha = function(data,estimate){
    alpha <- 1/var(data/estimate)
    alphaold = alpha + 5
    for(i in 1:10){
      alpha = alpha - neglikdiff1(alpha, data = data, estimate = estimate)/neglikdiff2(alpha, data = data, estimate = estimate)
      if(!(alpha > 0)){ alpha = runif(1,1,10)}

      if(abs(alpha - alphaold) < 0.01) break
      alphaold = alpha
    }
    return(alpha)
  }

  if(patient_specific){
    data <- t(data)
    alpha = numeric(ncol(data))
    estimate = t(w_p%*%h_p)
    for(i in 1:ncol(data)){
      alpha[i] = NR_alpha(data[,i],estimate[,i])
    }
  }else{
    data = as.vector(data)
    estimate = as.vector(w_p%*%h_p)

    alpha = NR_alpha(data,estimate)
  }

  return(alpha)
}
