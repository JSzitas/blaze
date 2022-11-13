#include <Rcpp.h>
using namespace Rcpp;
using namespace RcppEigen;

#include "arima.h"
#include "utils.h"

#include <RcppEigen.h>
// #include <Eigen/Dense>
// #include <Eigen/Map>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
void try_make_arima( std::vector<double> phi,
                     std::vector<double> theta,
                     std::vector<double> delta) {

  auto res = make_arima(phi, theta, delta);
  // checked:
  // print_vector(res.V);
  // print_vector(res.Z);
  // print_vector(res.a);
  // print_vector(res.P);
  // print_vec_as_sqr_mat(res.T);
  // TODO:
  // print_vec_as_sqr_mat(res.Pn);

  return;
}
