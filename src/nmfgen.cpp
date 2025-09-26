#include <RcppArmadillo.h>
#include <cmath>        
#include <tuple>
#include <iostream>
#include "Dev_functions.h"
#include "Update_functions.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List nmfLh(arma::mat data, arma::mat wmat, arma::mat hmat, int iter){
  arma::mat estimate;
  arma::mat absdiff;
  
  for (int i = 0; i < iter; i++){
    estimate = wmat * hmat;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    absdiff = arma::abs(data-estimate);
    absdiff.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    wmat = wmat % (((data/absdiff) * arma::trans(hmat))/((estimate/absdiff) * arma::trans(hmat)));
    wmat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    estimate = wmat * hmat;
    
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    absdiff = arma::abs(data-estimate);
    absdiff.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    hmat = hmat % ((arma::trans(wmat) * (data/absdiff))/(arma::trans(wmat) * (estimate/absdiff)));
    hmat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    
    }
  
  List output = List::create(Named("W") = wmat,
                             Named("H") = hmat);
  return output;
}

// [[Rcpp::export]]
List nmfall(arma::mat data, int noSignatures, std::string distribution, std::string method, arma::colvec alpha, double pwr, arma::mat wmat, arma::mat hmat, arma::mat w1mat, arma::mat w2mat, int maxiter = 10000, double tolerance = 1e-8, int initial = 100, int smallIter = 500) {
  std::tuple<arma::mat, arma::mat, double> res;
  arma::mat wmat0 = Rcpp::NumericMatrix::create();
  arma::mat hmat0 = Rcpp::NumericMatrix::create();
  arma::mat w1mat0 = Rcpp::NumericMatrix::create();
  arma::mat w2mat0 = Rcpp::NumericMatrix::create();
  double errorValue = R_NaN;
  
  if (distribution == "NegativeBinomial"){
    res = NBupdates(data, noSignatures, alpha, method, smallIter, 0, wmat, hmat, w1mat, w2mat);
  } else if (distribution == "Tweedie") {
    res = TWupdates(data, noSignatures, pwr,  method, smallIter, 0, wmat, hmat, w1mat, w2mat);
  } else if (distribution == "Laplace"){
    res = Lupdates(data, noSignatures, method, smallIter, 0, wmat, hmat, w1mat, w2mat);
  }

  if (method == "traditional"){
    wmat0 = std::get<0>(res);
    hmat0 = std::get<1>(res);
    errorValue = std::get<2>(res);
  } else {
    w1mat0 = std::get<0>(res);
    w2mat0 = std::get<1>(res);
    errorValue = std::get<2>(res);
  }

  for(int i = 1; i < initial; i++){
    if (distribution == "NegativeBinomial"){
      res = NBupdates(data, noSignatures, alpha, method, smallIter, 0, wmat, hmat, w1mat, w2mat);
    } else if (distribution == "Tweedie") {
      res = TWupdates(data, noSignatures, pwr,  method, smallIter, 0, wmat, hmat, w1mat, w2mat);
    } else if (distribution == "Laplace"){
      res = Lupdates(data, noSignatures, method, smallIter, 0, wmat, hmat, w1mat, w2mat);
    }
    auto errorNew = std::get<2>(res);
    
    if(errorNew < errorValue){
      errorValue = errorNew;
      if (method == "traditional"){
        wmat0 = std::get<0>(res);
        hmat0 = std::get<1>(res);
        w1mat0 = Rcpp::NumericMatrix::create();
        w2mat0 = Rcpp::NumericMatrix::create();
        errorValue = std::get<2>(res);
      } else {
        w1mat0 = std::get<0>(res);
        w2mat0 = std::get<1>(res);
        wmat0 = Rcpp::NumericMatrix::create();
        hmat0 = Rcpp::NumericMatrix::create();
        errorValue = std::get<2>(res);
      }
    }
  }

  if (distribution == "NegativeBinomial"){
    res = NBupdates(data, noSignatures, alpha, method, maxiter, tolerance, wmat0, hmat0, w1mat0, w2mat0);
  } else if (distribution == "Tweedie") {
    res = TWupdates(data, noSignatures, pwr,  method, maxiter, tolerance, wmat0, hmat0, w1mat0, w2mat0);
  } else if (distribution == "Laplace"){
    res = Lupdates(data, noSignatures, method, maxiter, tolerance, wmat0, hmat0, w1mat0, w2mat0);
  }
  
  if (method == "traditional"){
    wmat = std::get<0>(res);
    hmat = std::get<1>(res);
    errorValue = std::get<2>(res);
    
    arma::colvec rsum = sum(hmat,1);
    wmat = wmat.each_row() % arma::trans(rsum);
    hmat = hmat.each_col() / rsum;
    
    List output = List::create(Named("W") = wmat,
                               Named("H") = hmat,
                               Named("cost") = errorValue);
    return output;
  } else {
    w1mat = std::get<0>(res);
    w2mat = std::get<1>(res);
    errorValue = std::get<2>(res);
    
    arma::mat hmat = arma::trans(w1mat) * data;
    arma::colvec rsum = sum(hmat,1);
    w2mat = w2mat.each_col() % rsum;
    w1mat = w1mat.each_row() / arma::trans(rsum);
    
    List output = List::create(Named("W1") = w1mat,
                               Named("W2") = w2mat,
                               Named("cost") = errorValue);
    return output;
  }
  
}