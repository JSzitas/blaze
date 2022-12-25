#include <Rcpp.h>
using namespace Rcpp;
#include "arima.h"
#include "string"

// [[Rcpp::plugins(cpp17)]]


class BlazeArima {
public:
  BlazeArima( std::vector<double> y,
              std::vector<int> orders_period = {1,0,1,0,0,0,1},
              std::vector<std::vector<double>> xreg = {{}},
              std::string ss_init = std::string("Gardner"),
              std::vector<bool> intercept_transform = {true, true},
              double kappa = 1000000) {

    arima_kind kind = arima_kind(orders_period);
    // the map is certainly nicer than if else statements and switches
    std::map<std::string, SSinit> init_map;
    init_map["Gardner"] = SSinit::Gardner;
    init_map["Rossignol"] = SSinit::Rossignol;
    this->model = Arima<double>( y, kind, xreg, intercept_transform[0],
                                intercept_transform[1],
                                init_map[ss_init],
                                CSS,
                                kappa);
  }
  void fit(){
    this->model.fit();
  }
  std::vector<double> get_coef() {
    return this->model.get_coef();
  }
  Rcpp::List forecast( int h = 10, std::vector<std::vector<double>> newxreg = {{}} ){
    forecast_result<double> result = this->model.forecast(h, newxreg );
    return List::create(Named("forecast") = result.forecast, Named("std.err.") = result.se);
  };
private:
  Arima<double> model;
};

RCPP_EXPOSED_CLASS_NODECL(BlazeArima)
RCPP_MODULE(BlazeArima) {
  // using namespace Rcpp;

  Rcpp::class_<BlazeArima>("BlazeArima")
  .constructor< std::vector<double>,
                std::vector<int>,
                std::vector<std::vector<double>>,
                std::string,
                std::vector<bool>,
                double>("basic contructor")
  .method("fit", &BlazeArima::fit, "fit arima model")
  .method("get_coef", &BlazeArima::get_coef, "get fitted arima coefficients")
  .method("forecast", &BlazeArima::forecast, "forecast from a fitted arima model");
}
