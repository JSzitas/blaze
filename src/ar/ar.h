#ifndef BLAZE_AR
#define BLAZE_AR

#include "utils/utils.h"
#include "utils/xreg.h"

#include "common/forecast_result.h"

#include "ar/ar_fitter.h"

template <typename scalar_t, typename Scaler = StandardScaler<scalar_t>> class AR {
private:
  // data
  std::vector<scalar_t> y;
  const size_t p;
  std::vector<std::vector<scalar_t>> xreg;
  bool intercept, drift;
  // estimated during fitting
  std::vector<scalar_t> coef, residuals, fitted_vals, reg_coef;
  scalar_t sigma2;
  std::vector<Scaler> scalers;
  scalar_t aic;
  bool fitted;
public:
  AR<scalar_t, Scaler>(){};
  AR<scalar_t, Scaler>(
      const std::vector<scalar_t> &y,
      const size_t p,
      const std::vector<std::vector<scalar_t>> &xreg = {{}},
      const bool intercept = true, const bool drift = true,
      const bool standardize = true) : y(y), p(p), xreg(xreg),
      intercept(intercept), drift(drift){
    this->scalers = std::vector<Scaler>( standardize * (1 + xreg.size()) );
    this->fitted = false;
    this->residuals = std::vector<scalar_t>(y.size(), 0.0);
    this->fitted_vals = std::vector<scalar_t>(y.size(), 0.0);
  };
  void fit() {
    // if we are asked to fit a model with p,q, P,Q == 0, we simply throw 
    if(p == 0) {
      std::cout << "You specified a model with no parameters (p equal" <<
        " to zero) - we will not fit a model in this case. Please" <<
          " use the .respecify() method and provide a valid 'p' parameter." <<
            std::endl;
      return;
    }
    // this should just proceed with fitting, not do things which can be done in
    // the constructor
    // if we have any scalers, fit them and apply scaling
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0] = Scaler(this->y);
      this->scalers[0].scale(this->y);
      size_t i = 1;
      for( auto & xreg_val:this->xreg ) {
        this->scalers[i] = Scaler(xreg_val);
        this->scalers[i].scale(xreg_val);
        i++;
      }
    }
    this->coef = std::vector<scalar_t>(this->xreg.size() + this->p + this->intercept + this->drift);
    this->sigma2 = fit_ar(coef, this->fitted_vals, this->residuals, 
                          this->y, this->xreg, this->p, this->intercept, 
                          this->drift);
    // we have to rescale the sigma to be on the same scale as original data
    auto sigma_2 = this->sigma2;
    if( scalers.size() > 0 ) {
      sigma_2 = scalers[0].rescale_val_w_mean(this->sigma2);
    }
    // 1.837877 is equal to log(2*pi) - log is not standard compliant in a
    // constexpr so we must expand the expression manually, sadly
    constexpr scalar_t one_p_log_twopi = 1.0 + 1.837877;
    this->aic = (this->y.size() - this->p) * (log(sigma_2) + one_p_log_twopi);
    // invert scaling - this is probably redundant, we do not care for estimated
    // coefficients too much (as the main goal of the package is forecasting)
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0].rescale(this->y);
      size_t i = 1;
      for( auto & xreg_val:this->xreg ) {
        scalar_t temp = this->scalers[i].rescale_val(coef[p+i-1]);
        this->coef[p+i-1] = temp;
        i++;
      }
      if( this->intercept ) {
        // the intercept also has to be rescaled to have any meaning
        scalar_t temp = scalers[0].rescale_val(this->coef.back());
        this->coef.back() = temp;
      }
      // rescale residuals to be on the original scale
      scalers[0].rescale_w_sd(this->residuals);
    }
    this->fitted = true;
  };
  void respecify(const size_t p) {
    this->p = p;
  }
  forecast_result<scalar_t> forecast(
      const size_t h = 10,
      std::vector<std::vector<scalar_t>> newxreg = {{}}) {
    // validate xreg length
    if (!this->fitted || newxreg.size() != this->xreg.size())
      return forecast_result<scalar_t>(0);
    // otherwise run forecast
    // TODO: rewrite forecast method
    auto res = kalman_forecast(h, this->model);
    // and if using xreg or intercepts, integrate those
    if( newxreg.size() > 0 || this->intercept ) {
      if( scalers.size() > 0 ) {
        size_t i = 1;
        for( auto & xreg_val:newxreg ) {
          scalers[i].scale(xreg_val);
          i++;
        }
      }
      auto xreg_adjusted = predict(h, this->reg_coef, this->intercept,
                                   this->drift, newxreg);
      for (size_t i = 0; i < h; i++) {
        res.forecast[i] += xreg_adjusted[i];
      }
    }
    // scale standard errors
    for( size_t i = 0; i < h; i++ ) {
      res.std_err[i] = sqrt(res.std_err[i] * this->sigma2);
    }
    if( scalers.size() > 0 ) {
      scalers[0].rescale(res.forecast);
      scalers[0].rescale_w_sd(res.std_err);
    }
    return res;
  };
  const std::vector<scalar_t> get_coef() const { return this->coef; }
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const bool is_fitted() const { return this->fitted; }
};

#endif
