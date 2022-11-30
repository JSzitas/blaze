#include <Rcpp.h>
using namespace Rcpp;

#include "utils/utils.h"
#include "xreg.h"
#include "arima_utils.h"

#include "string"

// [[Rcpp::export]]
std::vector<double> lm_cpp( std::vector<double> &y,
                            std::vector<std::vector<double>> & xreg,
                            bool use_intercept = true) {

  std::vector<double> xreg_;
  if( xreg.size() >0 ) {
    xreg_ = std::vector<double>( xreg.size() * xreg[0].size() );
    int p = 0;
    for( int j = 0; j < xreg.size(); j++ ) {
      for( int i = 0; i < xreg[0].size(); i++) {
        xreg_[p] = xreg[j][i];
        p++;
      }
    }
  }
  else {
    xreg_ = std::vector<double>(0);
  }

  // return xreg_;

  auto res = xreg_coef( y, xreg_, use_intercept);
  return res.coef;
}

// [[Rcpp::export]]
int can_you_string( std::string input ) {
  std::cout << input << std::endl;
  return 0;
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
