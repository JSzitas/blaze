#include <Rcpp.h>
using namespace Rcpp;
// using namespace RcppEigen;

// #include "arima.h"
// #include "utils.h"
#include "xreg.h"

// #include <RcppEigen.h>
// #include <Eigen/Dense>
// #include <Eigen/Map>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
std::vector<double> lm_cpp( std::vector<double> &y,
                            std::vector<std::vector<double>> & xreg,
                            bool use_intercept = true) {

  std::vector<double> xreg_( xreg.size() * xreg[0].size() );
  int p = 0;
  for( int j = 0; j < xreg.size(); j++ ) {
    for( int i = 0; i < xreg[0].size(); i++) {
      xreg_[p] = xreg[j][i];
      p++;
    }
  }

  auto res = xreg_coef( y, xreg_, use_intercept);
  return res;
}




// void count_non_na( std::vector<double> y) {
//
//   int n = 0;
//   for( auto & el:y) {
//     if( !std::isnan(el)) {
//       n++;
//     }
//   }
//
//   std::cout << n;
//
//   return;
// }
