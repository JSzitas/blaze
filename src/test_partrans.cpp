#include <Rcpp.h>
using namespace Rcpp;
#include "arima_utils2.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
void partrans_test( std::vector<double> x) {

  std::cout << "Initial: " << std::endl;
  print_vector(x);

  std::vector<double> y = x;
  auto res = parameter_transform(y);
  std::cout << "After partrans 1: "<< std::endl;
  print_vector(res);

  parameter_transform2(y, 0, x.size());
  std::cout << "After partrans 2: "<< std::endl;
  print_vector(y);
}

// [[Rcpp::export]]
std::vector<double> arima_partrans_test( std::vector<double> x, std::vector<int> sarima_structure = {3,1,2,2,1,1,3};) {

  arima_kind kind(sarima_structure);

  arima_transform_parameters2(x, kind, true);
  return x;
}





