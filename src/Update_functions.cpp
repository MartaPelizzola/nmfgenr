#include <RcppArmadillo.h>
#include <cmath>        // std::abs
#include <tuple>
#include <iostream>
#include <Rcpp.h>
#include "Dev_functions.h"
#include "Update_functions.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double gammapower(double pwr){
  double gp;
  if (pwr > 1) {
    gp = 1/(pwr);
  }
  if ((pwr >= 0) && (pwr <= 1)) {
    gp = 1;
  }
  if (pwr < 0) {
    gp = 1/(1-pwr);
  }
  return gp;
}


std::tuple<arma::mat, arma::mat, double> TWupdates(arma::mat data, int noSignatures, double pwr, std::string method, int iter, double tol, arma::mat wmat, arma::mat hmat, arma::mat w1mat, arma::mat w2mat) {
  int datapts = data.n_rows; // N
  int datadims = data.n_cols; // M

  if (method == "traditional"){
    if (wmat.is_empty() || hmat.is_empty()){
      hmat = arma::mat(noSignatures, datadims, arma::fill::randu);
      wmat = arma::mat(datapts, noSignatures, arma::fill::randu);
    }

    arma::mat estimate = wmat * hmat;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    double error = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);
    double errorOld = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);
    double errorNew = 2*errorOld;

    for(int t = 0; t < iter; t++){
      wmat = wmat % pow(((data/pow(estimate, pwr)) * arma::trans(hmat))/(pow(estimate, 1-pwr) * arma::trans(hmat)), gammapower(pwr));
      wmat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = wmat * hmat;
      estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      hmat = hmat % pow((arma::trans(wmat) * (data/pow(estimate, pwr)))/(arma::trans(wmat) * pow(estimate, 1-pwr)), gammapower(pwr));
      hmat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = wmat * hmat;

      if (tol != 0){
        if(t - floor(t/10)*10 == 0){
          errorNew = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);

          if (2*std::abs(errorOld - errorNew)/(0.1 + std::abs(2*errorNew)) < tol){
            break;
          }
          errorOld = errorNew;
        }
      }
    }
    error = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);

    return {wmat, hmat, error};
  } else if (method == "convex") {
    if (w1mat.is_empty() || w2mat.is_empty()){
      w1mat = arma::mat(datapts, noSignatures, arma::fill::randu);
      w2mat = arma::mat(noSignatures, datapts, arma::fill::randu);
    }

    arma::mat estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    double error = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);
    double errorOld = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);
    double errorNew = 2*errorOld;

    for(int t = 0; t < iter; t++){
      w2mat = w2mat % pow((arma::trans(w1mat)*data * arma::trans(pow(estimate, -pwr) % data))/(arma::trans(w1mat)*data * arma::trans(pow(estimate, 1-pwr))),gammapower(pwr));

      estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
      estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      w1mat = w1mat % pow((data * arma::trans(pow(estimate, -pwr) % data) * arma::trans(w2mat))/(data * arma::trans(pow(estimate, 1-pwr)) * arma::trans(w2mat)),gammapower(pwr));

      w1mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
      w2mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;

      if (tol != 0){
        if(t - floor(t/10)*10 == 0){
          errorNew = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);

          if (std::abs(errorOld - errorNew)/(std::abs(2*errorNew)) < tol){
            break;
          }
          errorOld = errorNew;
        }
      }
    }
    error = unitdev(arma::vectorise(data),arma::vectorise(estimate),pwr);

    return {w1mat, w2mat, error};
  } else {
    Rcout << "method needs to be either traditional or convex";
    return {wmat, hmat, R_NaN};
  }
}

std::tuple<arma::mat, arma::mat, double> NBupdates(arma::mat data, int noSignatures, arma::colvec alpha, std::string method, int iter, double tol, arma::mat wmat, arma::mat hmat, arma::mat w1mat, arma::mat w2mat) {
  int datapts = data.n_rows; // N
  int datadims = data.n_cols; // M

  arma::mat alphamat(datapts, datadims, arma::fill::zeros);
  alphamat.each_col() = alpha;

  if (method == "traditional"){
    if (wmat.is_empty() || hmat.is_empty()){
      hmat = arma::mat(noSignatures, datadims, arma::fill::randu);
      wmat = arma::mat(datapts, noSignatures, arma::fill::randu);
    }

    arma::mat estimate = wmat * hmat;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    double error = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
    double errorOld = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
    double errorNew = 2*errorOld;

    for(int t = 0; t < iter; t++){
      wmat = wmat % (((data/estimate) * arma::trans(hmat))/(((alphamat + data)/(alphamat + estimate))* arma::trans(hmat)));

      estimate = wmat * hmat;
      estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      hmat = hmat % ((arma::trans(wmat) * (data/estimate))/(arma::trans(wmat)*((alphamat + data)/(alphamat + estimate))));

      wmat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
      hmat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = wmat * hmat;

      if (tol != 0){
        if(t - floor(t/10)*10 == 0){
          errorNew = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));

          if (2*std::abs(errorOld - errorNew)/(0.1 + std::abs(2*errorNew)) < tol){
            break;
          }
          errorOld = errorNew;
        }
      }
    }
    error = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));

    return {wmat, hmat, error};

  } else if (method == "convex") {
    if (w1mat.is_empty() || w2mat.is_empty()){
      w1mat = arma::mat(datapts, noSignatures, arma::fill::randu);
      w2mat = arma::mat(noSignatures, datapts, arma::fill::randu);
    }

    arma::mat estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    double error = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
    double errorOld = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));
    double errorNew = 2*errorOld;

    for(int t = 0; t < iter; t++){
      w2mat = w2mat % ((arma::trans(w1mat)*data * arma::trans(data/estimate)) / (arma::trans(w1mat)*data * arma::trans((data + alphamat)/(estimate + alphamat))));

      estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
      estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      w1mat = w1mat % ((data * arma::trans(data/estimate) * arma::trans(w2mat)) / (data * arma::trans((data+alphamat)/(estimate+alphamat)) * arma::trans(w2mat)));

      w1mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
      w2mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;

      if (tol != 0){
        if(t - floor(t/10)*10 == 0){
          errorNew = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));

          if (std::abs(errorOld - errorNew)/(std::abs(2*errorNew)) < tol){
            break;
          }
          errorOld = errorNew;
        }
      }
    }
    error = nbdev(arma::vectorise(data),arma::vectorise(estimate),arma::vectorise(alphamat));

    return {w1mat, w2mat, error};
  } else {
    Rcout << "method needs to be either traditional or convex";
    return {wmat, hmat, R_NaN};
  }
}

std::tuple<arma::mat, arma::mat, double> Lupdates(arma::mat data, int noSignatures, std::string method, int iter, double tol, arma::mat wmat, arma::mat hmat, arma::mat w1mat, arma::mat w2mat) {
  int datapts = data.n_rows; // N
  int datadims = data.n_cols; // M

  if (method == "traditional"){
    if (wmat.is_empty() || hmat.is_empty()){
      hmat = arma::mat(noSignatures, datadims, arma::fill::randu);
      wmat = arma::mat(datapts, noSignatures, arma::fill::randu);
    }

    arma::mat estimate = wmat * hmat;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
    arma::mat absdiff = arma::abs(data-estimate);
    absdiff.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );


    double error = l1dev(arma::vectorise(data),arma::vectorise(estimate));
    double errorOld = l1dev(arma::vectorise(data),arma::vectorise(estimate));
    double errorNew = 2*errorOld;

    for(int t = 0; t < iter; t++){
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


      if (tol != 0){
        if(t - floor(t/10)*10 == 0){
          errorNew = l1dev(arma::vectorise(data),arma::vectorise(estimate));

          if (2*std::abs(errorOld - errorNew)/(0.1 + std::abs(2*errorNew)) < tol){
            break;
          }
          errorOld = errorNew;
        }
      }
    }
    error = l1dev(arma::vectorise(data),arma::vectorise(estimate));

    return {wmat, hmat, error};

  } else if (method == "convex") {
    if (w1mat.is_empty() || w2mat.is_empty()){
      w1mat = arma::mat(datapts, noSignatures, arma::fill::randu);
      w2mat = arma::mat(noSignatures, datapts, arma::fill::randu);
    }
    w2mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    arma::mat estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
    estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    arma::mat absdiff = arma::abs(data-estimate);
    absdiff.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

    double error = l1dev(arma::vectorise(data),arma::vectorise(estimate));
    double errorOld = l1dev(arma::vectorise(data),arma::vectorise(estimate));
    double errorNew = 2*errorOld;

    for(int t = 0; t < iter; t++){
      w1mat = w1mat % ((data * arma::trans(data/absdiff) * arma::trans(w2mat))/(data * arma::trans(estimate/absdiff) * arma::trans(w2mat)));
      w1mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
      estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );
      absdiff = arma::abs(data-estimate);
      absdiff.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      w2mat = w2mat % ((arma::trans(w1mat) * data * arma::trans(data/absdiff))/(arma::trans(w1mat) * data * arma::trans(estimate/absdiff)));
      w2mat.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      estimate = arma::trans(w2mat) * arma::trans(w1mat) * data;
      estimate.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );

      absdiff = arma::abs(data-estimate);
      absdiff.transform( [](double val) {return (val < 1e-16) ? 1e-16 : val; } );


      if (tol != 0){
        if(t - floor(t/10)*10 == 0){
          errorNew = l1dev(arma::vectorise(data),arma::vectorise(estimate));

          if (std::abs(errorOld - errorNew)/(std::abs(2*errorNew)) < tol){
            break;
          }
          errorOld = errorNew;
        }
      }
    }
    error = l1dev(arma::vectorise(data),arma::vectorise(estimate));

    return {w1mat, w2mat, error};
  } else {
    Rcout << "method needs to be either traditional or convex";
    return {wmat, hmat, R_NaN};
  }
}
