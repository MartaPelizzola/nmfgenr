#include <RcppArmadillo.h>
#include <cmath>        
#include <tuple>
#include <iostream>
using namespace Rcpp;
#include "Dev_functions.h"

// [[Rcpp::depends(RcppArmadillo)]]

double gkldev(arma::colvec y, arma::colvec mu) {
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
    if (y[i] <= 0 || mu[i] <= 0) {
      sum += mu[i];
    } else {
      sum += y[i] * (log(y[i]) - log(mu[i])) - y[i] + mu[i];
    }
  }
  return sum;
}

double unitdev(arma::colvec y, arma::colvec mu, double pwr){
  int ySize = y.size();
  double sum = 0;
  
  for (int i=0; i<ySize; i++) {
    if (y[i] > 0) {
      if (pwr==2) {
        sum += y[i]/mu[i] - log(y[i]) + log(mu[i]) - 1;
      } else if (pwr == 1) {
        sum += y[i] * (log(y[i]) - log(mu[i])) - y[i] + mu[i];
      } else {
        sum += ((pow(y[i], 2-pwr) + (1 - pwr)*pow(mu[i], 2-pwr) - (2-pwr)*y[i]*pow(mu[i], 1 -pwr))/((2-pwr)*(1-pwr)));
      }
    }
  }
  return sum;
}

double nbdev(arma::colvec y, arma::colvec mu, arma::colvec alpha){
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
    if (y[i] <= 0 || mu[i] <= 0) {
      sum += -(y[i] + alpha[i]) * (log(alpha[i] + y[i]) - log(alpha[i] + mu[i]));
    } else {
      sum += y[i] * (log(y[i]) - log(mu[i])) - (alpha[i] + y[i]) * (log(alpha[i] + y[i]) - log(alpha[i] + mu[i]));
    }
  }
  return sum;
}

double l1dev(arma::colvec y, arma::colvec mu){
  int ySize = y.size();
  double sum = 0;
  for (int i=0; i<ySize; i++) {
    sum += std::abs(y[i] - mu[i]);
  }
  return sum;
}