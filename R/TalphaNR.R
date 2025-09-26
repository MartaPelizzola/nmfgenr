#' @importFrom stats var
#' @importFrom stats runif
#'
#' @title Estimation of the overdispersion parameter for Negative Binomial non-negative matrix factorization from the SigMoS package: https://github.com/MartaPelizzola/SigMoS/tree/master
#'
#' @description Likehood estimation of the dispersion parameter in negative binomial using Newton-Raphson.
#' The overdispersion parameter can be either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#'
#'
#' @param data Numeric matrix of mutational counts data. Matrix size: no. of patients x no. of mutation types.
#' @param k Number of signatures to be used for the non-negative matrix factorization
#' @param patient_specific Logical. If TRUE patient-specific overdispersion is used in the Negative Binomial model.
#'
#'
#' @return Overdispersion parameter. Either vector of length one or vector of length no. of patients if patient-specific overdispersion is used.
#'
#' @examples
#' # Estimate patient specific overdispersion parameters:
#' alpha <- alphaNR(BRCA21,k=3, patient_specific = TRUE)
#'
#' @export
#'
#'
alphaNR <- function(data, k=NULL, patient_specific = FALSE){
}
