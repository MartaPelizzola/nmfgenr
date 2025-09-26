#ifndef UPDATE_FUNCTIONS_H
#define UPDATE_FUNCTIONS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <tuple>

using namespace Rcpp;

double gammapower(double pwr);
std::tuple<arma::mat, arma::mat, double> TWupdates(arma::mat data, int noSignatures, double pwr, std::string method, int iter, double tol = 0, arma::mat wmat = Rcpp::NumericMatrix::create(), arma::mat hmat = Rcpp::NumericMatrix::create(), arma::mat w1mat = Rcpp::NumericMatrix::create(), arma::mat w2mat = Rcpp::NumericMatrix::create()); 
std::tuple<arma::mat, arma::mat, double> NBupdates(arma::mat data, int noSignatures, arma::colvec alpha, std::string method, int iter, double tol = 0, arma::mat wmat = Rcpp::NumericMatrix::create(), arma::mat hmat = Rcpp::NumericMatrix::create(), arma::mat w1mat = Rcpp::NumericMatrix::create(), arma::mat w2mat = Rcpp::NumericMatrix::create()); 
std::tuple<arma::mat, arma::mat, double> Lupdates(arma::mat data, int noSignatures, std::string method, int iter, double tol = 0, arma::mat wmat = Rcpp::NumericMatrix::create(), arma::mat hmat = Rcpp::NumericMatrix::create(), arma::mat w1mat = Rcpp::NumericMatrix::create(), arma::mat w2mat = Rcpp::NumericMatrix::create()); 

#endif