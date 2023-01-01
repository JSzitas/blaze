#include <Rcpp.h>
using namespace Rcpp;

#include "arima/solvers/initializers.h"

// [[Rcpp::export]]
std::vector<double> test_get_Q0( std::vector<double> phi, std::vector<double> theta) {
  return get_Q0(phi, theta);
}


