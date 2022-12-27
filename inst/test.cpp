#include <Rcpp.h>
using namespace Rcpp;
#include "arima/solvers/state_space.h"
#include "arima/structures/arima_kind.h"

// [[Rcpp::plugins(cpp17)]]


// [[Rcpp::export]]
Rcpp::List test_make_arima( std::vector<double> phi,
                            std::vector<double> theta,
                            std::vector<double> delta) {//,
                    // std::vector<int> orders_period = {1,0,1,0,0,0,1}) {

  // arima_kind kind = arima_kind(orders_period);

  structural_model<double> res = make_arima(phi, theta, delta);
  return List::create(Named("phi") = res.phi,
                      Named("theta") = res.theta,
                      Named("Delta") = res.delta,
                      Named("Z") = res.Z,
                      Named("a") = res.a,
                      Named("P") = res.P,
                      Named("T") = res.T,
                      Named("V") = res.V,
                      Named("h") = res.h,
                      Named("Pn") = res.Pn);
}


// [[Rcpp::export]]
Rcpp::List kalman_forecast_test( int h_ahead, Rcpp::List model ) {//,

  std::vector<double> phi = model["phi"];
  std::vector<double> theta = model["theta"];
  std::vector<double> delta = model["Delta"];
  std::vector<double> Z = model["Z"];
  std::vector<double> a = model["a"];
  std::vector<double> P = model["P"];
  std::vector<double> T = model["T"];
  std::vector<double> V = model["V"];
  std::vector<double> Pn = model["Pn"];
  double h = model["h"];

  structural_model<double> mdl( phi, theta, delta, Z, a, P, T, V, h, Pn );
  forecast_result<double> result = kalman_forecast(h_ahead, mdl );
  return List::create(Named("forecast") = result.forecast, Named("std.err.") = result.std_err);
}
