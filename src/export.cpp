#include <Rcpp.h>
using namespace Rcpp;
#include "arima.h"
#include "string"

// [[Rcpp::plugins(cpp17)]]


class BlazeArima {
public:
  BlazeArima( std::vector<double> y,
              std::vector<int> orders_period,
              std::vector<std::vector<double>> xreg,
              std::string ss_init,
              std::vector<bool> intercept_transform,
              double kappa) {

    arima_kind kind = arima_kind(orders_period[0], orders_period[1], orders_period[2],
                                 orders_period[3], orders_period[4], orders_period[5],
                                 orders_period[6]);
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
    model.fit();
  }
  // Rcpp::List predict(){};
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
  .method("fit", &BlazeArima::fit, "fit arima model");
}
