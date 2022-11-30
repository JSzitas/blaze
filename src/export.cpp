#include <Rcpp.h>
using namespace Rcpp;
#include "arima.h"
#include "string"

class BlazeArima {
public:
  BlazeArima( std::vector<float> y,
              std::vector<int> orders_period,
              std::vector<std::vector<float>> xreg,
              std::string ss_init,
              std::vector<bool> intercept_transform,
              float kappa) {

    arima_kind kind = arima_kind(orders_period[0], orders_period[1], orders_period[2],
                                 orders_period[3], orders_period[4], orders_period[5],
                                 orders_period[6]);
    // the map is certainly nicer than if else statements and switches
    std::map<std::string, SSinit> init_map;
    init_map["Gardner"] = SSinit::Gardner;
    init_map["Rossignol"] = SSinit::Rossignol;
    this->model = Arima<float>( y, kind, xreg, intercept_transform[0],
                                intercept_transform[1], init_map[ss_init],
                                kappa);
  }
  void fit(){}
  // Rcpp::List predict(){};
private:
  Arima<float> model;
};

RCPP_EXPOSED_CLASS_NODECL(BlazeArima)
RCPP_MODULE(BlazeArima) {
  // using namespace Rcpp;

  Rcpp::class_<BlazeArima>("BlazeArima")
  .constructor< std::vector<float>,
                std::vector<int>,
                std::vector<std::vector<float>>,
                std::string,
                std::vector<bool>,
                float>("basic contructor");
  // .method("fit", &arima::fit, "fit arima model");
}
