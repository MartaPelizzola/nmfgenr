#ifndef DEV_FUNCTIONS_H
#define DEV_FUNCTIONS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <tuple>
using namespace Rcpp;

double gkldev(arma::colvec y, arma::colvec mu);
double unitdev(arma::colvec y, arma::colvec mu, double pwr);
double nbdev(arma::colvec y, arma::colvec mu, arma::colvec alpha);
double l1dev(arma::colvec y, arma::colvec mu);

#endif