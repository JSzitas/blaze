#include <Rcpp.h>
using namespace Rcpp;
#include "arima/arima.h"
#include "ar/ar.h"
#include "auto_ar/auto_ar.h"
#include "inoise/inoise.h"
#include "string"

// [[Rcpp::plugins(cpp17)]]

class arima {
private:
  Arima<double, StandardScaler<double>> model;
public:
  arima(std::vector<double> y,
        std::vector<size_t> orders_period = {1,0,1,0,0,0,1},
        std::vector<std::vector<double>> xreg = {{}},
        std::string ss_init = std::string("Gardner"),
        std::string method = std::string("CSS"),
        std::vector<bool> intercept_drift_transform = {true, false, true},
        double kappa = 1000000) {

    const arima_kind kind = arima_kind(orders_period);
    // map is certainly nicer than if else statements and switches
    static std::map<std::string, SSinit> ss_method_map{
      {"Gardner", SSinit::Gardner},
      {"Rossignol", SSinit::Rossignol}
    };
    static std::map<std::string, fitting_method> fitting_method_map = {
      {"CSS", fitting_method::CSS},
      {"CSS-ML", fitting_method::CSSML},
      {"ML", fitting_method::ML},
      {"LASSO", fitting_method::LASSO}
    };

    this->model = Arima<double, StandardScaler<double>>(
      y, kind, xreg, intercept_drift_transform[0],
      intercept_drift_transform[1], intercept_drift_transform[2],
      ss_method_map[ss_init], fitting_method_map[method], kappa, true);
  }
  void fit() {
    this->model.fit();
  }
  std::vector<double> get_coef() const {
    return this->model.get_coef();
  }
  std::vector<double> residuals() const {
    return this->model.get_residuals();
  };
  std::vector<double> fitted() const {
    return this->model.get_fitted();
  };
  // R specific exports 
  Rcpp::List get_structural_model() const {
    auto res = this->model.get_structural_model();
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
  Rcpp::List forecast( size_t h = 10, std::vector<std::vector<double>> newxreg = {{}} ){
    forecast_result<double> result = this->model.forecast(h, newxreg );
    return List::create(Named("forecast") = result.forecast, Named("std.err.") = result.std_err);
  };
};

class ar {
private:
  AR<double, StandardScaler<double>> model;
public:
  ar(const std::vector<double> &y,
     const size_t p,
     const std::vector<std::vector<double>> &xreg = {{}},
     const bool intercept = true, const bool drift = true) {

    this->model = AR<double, StandardScaler<double>>(
      y, p, xreg, intercept, drift, true);
  }
  void fit() {
    this->model.fit();
  }
  std::vector<double> get_coef() const {
    return this->model.get_coef();
  }
  std::vector<double> residuals() const {
    return this->model.get_residuals();
  };
  std::vector<double> fitted() const {
    return this->model.get_fitted();
  };
  // R specific exports 
  Rcpp::List forecast( size_t h = 10, std::vector<std::vector<double>> newxreg = {{}} ){
    forecast_result<double> result = this->model.forecast(h, newxreg );
    return List::create(Named("forecast") = result.forecast, Named("std.err.") = result.std_err);
  };
};

class auto_ar {
private:
  AutoAR<double, StandardScaler<double>> model;
public:
  auto_ar(const std::vector<double> &y,
     const size_t min_p,
     const size_t max_p,
     const std::vector<std::vector<double>> &xreg = {{}},
     const bool intercept = true, const bool drift = true, 
     std::string method = std::string("AIC")) {

    static std::map<std::string, AutoARMethod> fitting_method_map = {
      {"AIC", AutoARMethod::AIC},
      {"AICc", AutoARMethod::AICc},
      {"BIC", AutoARMethod::BIC}
    };

    this->model = AutoAR<double, StandardScaler<double>>(
      y, min_p, max_p, xreg, intercept, drift, true, fitting_method_map[method]);
  }
  void fit() {
    this->model.fit();
  }
  std::vector<double> get_coef() const {
    return this->model.get_coef();
  }
  std::vector<double> residuals() const {
    return this->model.get_residuals();
  };
  std::vector<double> fitted() const {
    return this->model.get_fitted();
  };
  // R specific exports 
  Rcpp::List forecast( size_t h = 10, std::vector<std::vector<double>> newxreg = {{}} ){
    forecast_result<double> result = this->model.forecast(h, newxreg);
    return List::create(Named("forecast") = result.forecast, Named("std.err.") = result.std_err);
  };
};

class inoise {
private:
  INoise<double, sqrtStabilizer<double>> model;
public:
  inoise() = delete;
  inoise(const std::vector<double> &y,
         const size_t differences = 1,
         const size_t seasonal_difference = 10,
         const bool stabilize = true) {
    this->model = INoise<double,sqrtStabilizer<double>>(
      y, differences, seasonal_difference, stabilize);
  }
  void fit() {
    this->model.fit();
  }
  // R specific exports 
  Rcpp::List forecast(const size_t h = 10, const bool produce_se = true,
                      const size_t samples = 10000 ){
    forecast_result<double> result = this->model.forecast(h, produce_se, samples );
    return List::create(Named("forecast") = result.forecast, Named("std.err.") = result.std_err);
  };
};


RCPP_EXPOSED_CLASS_NODECL(arima)
RCPP_MODULE(blaze_arima) {
  Rcpp::class_<arima>("blaze_arima")
  .constructor< std::vector<double>,
                std::vector<size_t>,
                std::vector<std::vector<double>>,
                std::string,
                std::string,
                std::vector<bool>,
                double>("contructor")
  .method("fit", &arima::fit, "fit arima model")
  .method("get_coef", &arima::get_coef, "get fitted arima coefficients")
  .method("forecast", &arima::forecast, "forecast from a fitted arima model")
  .method("get_structural_model", &arima::get_structural_model,
                "get structural model specification from arima model")
  .method("residuals", &arima::residuals, "get residuals from fitted model")
  .method("fitted", &arima::fitted, "get fitted values from fitted model");
}

RCPP_EXPOSED_CLASS_NODECL(ar)
RCPP_MODULE(blaze_ar) {
  Rcpp::class_<ar>("blaze_ar")
  .constructor<std::vector<double>, size_t, std::vector<std::vector<double>>,
               bool, bool>("contructor")
  .method("fit", &ar::fit, "fit ar model")
  .method("get_coef", &ar::get_coef, "get fitted ar coefficients")
  .method("forecast", &ar::forecast, "forecast from a fitted ar model")
  .method("residuals", &ar::residuals, "get residuals from fitted model")
  .method("fitted", &ar::fitted, "get fitted values from fitted model");
}

RCPP_EXPOSED_CLASS_NODECL(auto_ar)
  RCPP_MODULE(blaze_auto_ar) {
    Rcpp::class_<auto_ar>("blaze_auto_ar")
    .constructor<std::vector<double>, size_t, size_t, std::vector<std::vector<double>>,
    bool, bool, std::string>("contructor")
    .method("fit", &auto_ar::fit, "fit ar model")
    .method("get_coef", &auto_ar::get_coef, "get fitted ar coefficients")
    .method("forecast", &auto_ar::forecast, "forecast from a fitted ar model")
    .method("residuals", &auto_ar::residuals, "get residuals from fitted model")
    .method("fitted", &auto_ar::fitted, "get fitted values from fitted model");
  }

RCPP_EXPOSED_CLASS_NODECL(inoise)
  RCPP_MODULE(blaze_inoise) {
    Rcpp::class_<inoise>("blaze_inoise")
    .constructor<std::vector<double>, size_t, size_t, bool>("contructor")
    .method("fit", &inoise::fit, "fit ar model")
    .method("forecast", &inoise::forecast, "forecast from a fitted ar model");
  }
