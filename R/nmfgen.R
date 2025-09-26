#' A general purpose function to perform traditional and convex NMF with different divergence measures.
#'
#' @param data A matrix of count data.
#' @param rank The rank for the factorization.
#' @param distribution The underlying model assumption for the factorization. Can assume the following values: "NegativeBinomial", "Tweedie", or "Laplace".
#' @param method The chosen method for the factorization. Can be either "traditional" or "convex".
#' @param alpha Dispersion parameter for the Negative Binomial matrix factorization.
#' @param pwr Power parameter for the Tweedie matrix factorization.
#' @param wmat An initialization value for the W matrix for traditional NMF.
#' @param hmat An initialization value for the H matrix for traditional NMF.
#' @param w1mat An initialization value for the W1 matrix for convex NMF.
#' @param w2mat An initialization value for the W2 matrix for convex NMF.
#' @param maxiter The maximum number of iterations for the estimation.
#' @param tolerance The tolerance value for convergence.
#' @param initial The number of initializations to test.
#' @param smallIter The number of iterations for which to run each initialization.
#'
#' @returns A list with the values for the estimated matrices (W and H for traditional NMF or W1 and W2 for convex NMF) and the final cost value.
#' @export
#'
#' @examples
nmfgen <- function(data, rank, distribution = c("NegativeBinomial", "Tweedie", "Laplace"), method = c("traditional", "convex"), alpha, pwr, wmat, hmat, w1mat, w2mat, maxiter, tolerance, initial, smallIter){
  if (is.null(data) | length(data)!=nrow(data)*ncol(data)){
    stop("data must be a matrix.")
  }
  if (is.null(rank) | !is.integer(rank) | length(rank)>1){
    stop("rank must be provided and needs to be a single integer number.")
  }
  if (method == "NegativeBinomial" & is.null(alpha)){
    alpha <- convex_alphaNR(data, rank, TRUE)
  }
  if (method == "Tweedie" & is.null(pwr)){
    pwr <- tweedie.profile()
  }
  if (!is.null(wmat) & (nrow(wmat)!=nrow(data) | ncol(wmat)!=rank)){
    stop("For traditional NMF the W matrix needs to have the same number of rows as the data matrix and 'rank' columns.")
  }
  if (!is.null(hmat) & (nrow(hmat)!=rank | ncol(hmat)!=ncol(data))){
    stop("For traditional NMF the H matrix needs to have 'rank' rows and the same number of columns as the data matrix.")
  }
  if (!is.null(w1mat) & (nrow(w1mat)!=nrow(data) | ncol(w1mat)!=rank)){
    stop("For traditional NMF the W1 matrix needs to have the same number of rows as the data matrix and 'rank' columns.")
  }
  if (!is.null(w2mat) & (nrow(w2mat)!=nrow(data) | ncol(w2mat)!=rank)){
    stop("For traditional NMF the W2 matrix needs to have the same number of rows as the data matrix and 'rank' columns.")
  }
  if (is.null(maxiter)){
    maxiter = 100000
    warning("maxiter not provided, the estimation will run for up to 10000 iterations.")
  }
  if (is.null(tolerance)){
    tolerance = 1e-8
    warning("tolerance not provided, convergence will be defined up to a tolerance of 1e-8.")
  }
  if (is.null(initial)){
    initial = 50
    warning("initial not provided, 100 initializations values will be tested.")
  }
  if (is.null(smallIter)){
    smallIter = 100
    warning("smallIter not provided, the initialization runs will run for 500 iterations each.")
  }
  if(is.null(wmat)){
    wmat = matrix(numeric(0))
  }
  if(is.null(hmat)){
    hmat = matrix(numeric(0))
  }
  if(is.null(w1mat)){
    w1mat = matrix(numeric(0))
  }
  if(is.null(w2mat)){
    w2mat = matrix(numeric(0))
  }
  if (distribution == "NegativeBinomial"){
    res <- nmfall(data, rank, distribution, method, alpha, NULL, wmat, hmat, w1mat, w2mat, maxiter, tolerance, initial, smallIter)
  }
  if (distribution == "Tweedie"){
    res <- nmfall(data, rank, distribution, method, NULL, pwr, wmat, hmat, w1mat, w2mat, maxiter, tolerance, initial, smallIter)
  }
  if (distribution == "Laplace"){
    res <- nmfall(data, rank, distribution, method, NULL, NULL, wmat, hmat, w1mat, w2mat, maxiter, tolerance, initial, smallIter)
  }
  return(res)
}
