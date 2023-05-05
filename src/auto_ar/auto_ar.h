#ifndef BLAZE_AUTO_AR
#define BLAZE_AUTO_AR

// #include "utils/utils.h"
// #include "utils/xreg.h"
// 
// #include "common/forecast_result.h"
// 
// #include "ar/ar_fitter.h"
#include "ar/ar.h"

enum AutoARMethod {
  AIC
};

template <typename scalar_t, typename Scaler = StandardScaler<scalar_t>> class AutoAR {
private:
  // data
  std::vector<scalar_t> y;
  size_t min_p, max_p;
  std::vector<std::vector<scalar_t>> xreg;
  bool intercept, drift, standardize;
  AutoARMethod method;
  // fitted best model
  AR<scalar_t, Scaler> fitted_ar;
public:
  AutoAR<scalar_t, Scaler>(){};
  AutoAR<scalar_t, Scaler>(
      const std::vector<scalar_t> &y,
      const size_t min_p,
      const size_t max_p,
      const std::vector<std::vector<scalar_t>> &xreg = {{}},
      const bool intercept = true, const bool drift = true,
      const bool standardize = true, 
      AutoARMethod method = AIC) : y(y), min_p(min_p), max_p(max_p),
      xreg(xreg), intercept(intercept), drift(drift), standardize(standardize),
      method(method){};
  void fit() {
    if(this->method == AutoARMethod::AIC) {
      this->fit_w_aic();
    }
    else {
      std::cout << "Invalid fitting method - please use respecify to provide " <<
        "a valid method, ideally one of " << 
        "AIC" << std::endl;
    }
  }
  void fit_w_aic() {
    // if we are asked to fit a model with p,q, P,Q == 0, we simply throw 
    if(max_p == 0 || min_p > max_p) {
      std::cout << "You specified a model with no parameters (max_p equal" <<
        " to zero) or a model where min_p > max_p - we will not fit a model " <<
          " in this case. Please use the .respecify() method and provide  " <<
            "valid 'min_p' and 'max_p' parameters." <<
            std::endl;
      return;
    }
    this->fitted_ar = AR<scalar_t, Scaler>(
      this->y, min_p, this->xreg, this->intercept, this->drift, this->standardize);
    this->fitted_ar.fit();
    // get aic from first model
    scalar_t best_aic = this->fitted_ar.get_aic(); 
    for( size_t p = this->min_p+1; p <= max_p; p++ ) {
      // create candidate model
      AR<scalar_t, Scaler> candidate(
          this->y, p, this->xreg, this->intercept, this->drift, this->standardize);
      // fit
      candidate.fit();
      scalar_t current_aic = candidate.get_aic();
      // check aic
      if(current_aic < best_aic) {
        // if best, set the fitted model to candidate 
        this->fitted_ar = candidate;
        // and update best aic
        best_aic = current_aic;
      }
    }
  };
  void respecify(const size_t min_p, const size_t max_p) {
    this->min_p = min_p;
    this->max_p = max_p;
  }
  void respecify(AutoARMethod method) {
    this->method = method;
  }
  forecast_result<scalar_t> forecast(
      const size_t h = 10,
      std::vector<std::vector<scalar_t>> &newxreg = {{}}) {
    return this->fitted_ar.forecast(h, newxreg);
  };
  const std::vector<scalar_t> get_coef() const { return this->fitted_ar.get_coef(); }
  const std::vector<scalar_t> get_residuals() const { return this->fitted_ar.get_residuals(); }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_ar.get_fitted(); }
  const scalar_t get_sigma2() const { return this->fitted_ar.get_sigma2(); }
  const scalar_t get_aic() const { return this->fitted_ar.get_aic(); }
  const bool is_fitted() const { return this->fitted_ar.fitted; }
};

#endif
