#' Profile likelihood estimation of power for the Tweedie distribution
#'
#' @param data A matrix of count data.
#' @param rank The rank for the factorization.
#' @param method The chosen method for the factorization. Can be either "traditional" or "convex".
#' @param power.vector A numeric vector with values between 1 and 2 or equal to 0.
#'
#' @returns A list with the estimated power parameter of the Tweedie distribution and the full results from the estimation procedure.
#' @export
#'
#' @examples

powertweedie <- function(data,  rank, method = c("traditional", "convex"), power.vector) {
  if(2 %in% power.vector && 0 %in% data){
    break("The gamma distribution is unsuitable for sparse data")
  }
  if(sum((power.vector < 0) | (0 > power.vector & power.vector < 1) | (power.vector > 2)) > 0){
    break("The values of power must be in the union of 0 and [1,2]")
  }

  x <- as.vector(data)
  phi <- c()
  loglik <- c()
  convergencecheck <- c()

  for(i in 1:length(power.vector)){
    if(method == "traditional"){
      nmf <- nmfgen(data, rank , "Tweedie", method = "traditional", pwr = power.vector[i])
      estimate <- nmf$W %*% nmf$H
    }else if(nmfmodel %in% c("Convex", "convex", "C")){
      nmf <- nmfgen(data, rank, "Tweedie", method = "convex", pwr = power.vector[i])
      estimate <- t(nmf$W2) %*% t(nmf$W1) %*% data
    }

    fit <- glm(x ~ as.vector(estimate), family = tweedie(var.power = power.vector[i], link.power = 1), start = c(0,1))


    phi[i] <- sigma(fit)^2

    if(power.vector[i] == 0){
      loglik[i] <- sum(dnorm(x, mean = as.vector(estimate), sd = sigma(fit), log = TRUE))
    }else if(power.vector[i] == 1){
      loglik[i] <- sum(ldTweedie(x, mu = as.vector(estimate), p = power.vector[i], phi = 1)[,1])
    }else{
      loglik[i] <- sum(ldTweedie(x, mu = as.vector(estimate), p = power.vector[i], phi = sigma(fit)^2)[,1])
    }
  }

  best <- which.max(loglik)
  fullresults <- list(log.likelihood = loglik,
                      phi.vec = phi,
                      power.vec = power.vector)

  results <- list(best.power = power.vector[best],
                  best.phi = phi[best],
                  best.loglik = loglik[best],
                  more = fullresults)
  return(results)
}
