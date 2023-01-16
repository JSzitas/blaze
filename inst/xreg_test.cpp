#include <Rcpp.h>
using namespace Rcpp;

#include "arima/utils/xreg.h"

// [[Rcpp::export]]
std::vector<double> test_xreg(std::vector<double> y,
                              std::vector<std::vector<double>> xreg,
                              const size_t d = 0,
                              const size_t D = 0,
                              const size_t period = 1,
                              bool intercept = true) {
  lm_coef<double> reg_coef( xreg.size(), intercept );;
  // fit xreg
  if (xreg.size() > 0 || intercept) {
    std::vector<double> y_d;
    std::vector<std::vector<double>> xreg_d;
    // if we have any differences
    if (d > 0) {
      y_d = diff(y, 1, d);
      xreg_d = diff(xreg, 1, d);
    }
    // seasonal differences
    if (period > 1 && D > 0) {
      y_d = diff(y, period, D);
      xreg_d = diff(xreg, period, D);
    }
    // fit coefficients and adjust y for fitted coefficients -
    // the original R code does this repeatedly, but it can be done only once
    // - the fitted effects from external regressors are never refitted
    if (y_d.size() <= xreg_d.size()) {
      reg_coef = xreg_coef(y, xreg, intercept);
    } else {
      reg_coef = xreg_coef(y_d, xreg_d);
    }
  }
  return reg_coef.data();
}

