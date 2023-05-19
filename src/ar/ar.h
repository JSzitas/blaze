#ifndef BLAZE_AR
#define BLAZE_AR

#include "utils/utils.h"
#include "utils/xreg.h"

#include "common/forecast_result.h"

#include "ar/ar_fitter.h"
#include "ar/ar_forecaster.h"

template <typename scalar_t, typename Scaler = StandardScaler<scalar_t>> class AR {
private:
  // data
  std::vector<scalar_t> y;
  size_t p;
  std::vector<std::vector<scalar_t>> xreg;
  bool intercept, drift;
  // estimated during fitting
  std::vector<scalar_t> coef, raw_coef, residuals, fitted_vals, prev_y;
  scalar_t sigma2;
  std::vector<Scaler> scalers;
  scalar_t aic, bic, aicc;
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
    this->prev_y = std::vector<scalar_t>(p+1);
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
    for(size_t i = 0; i < this->prev_y.size(); i++) {
      const size_t index = this->y.size() - this->prev_y.size() + i;
      this->prev_y[i] = this->y[index];
    }
    this->coef = std::vector<scalar_t>(
      this->xreg.size() + this->p + this->intercept + this->drift);
    this->raw_coef = std::vector<scalar_t>(this->coef.size());
    this->sigma2 = fit_ar(coef, this->fitted_vals, this->residuals, 
                          this->y, this->xreg, this->p, this->intercept, 
                          this->drift);
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0].rescale(this->y);
      size_t i = 1;
      // write back raw coefficients - we will need them later
      this->raw_coef = this->coef;
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
      // rescale fitted values 
      this->scalers[0].rescale(this->fitted_vals);
      for(size_t i = 0; i < p; i++) {
        // first p values set to zero
        this->fitted_vals[i] = 0.0;
      }
      // rescale residuals to be on the original scale
      this->scalers[0].rescale_w_sd(this->residuals);
    }
    this->aic = this->y.size() * std::log(
      crossprod(this->residuals)/this->y.size()) + 
      (2 * this->coef.size());
    // these need to be verified (they might be slightly off)
    this->bic = this->aic + (this->p + 1) * (log(this->y.size()) - 2);
    this->aicc = this->aic + 
      (2*(std::pow(this->p, 2) + this->p)/(this->y.size() - this->p - 1)); 
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
    // scale xreg as needed
    if( newxreg.size() > 0 || this->intercept ) {
      if( scalers.size() > 0 ) {
        size_t i = 1;
        for( auto & xreg_val:newxreg ) {
          scalers[i].scale(xreg_val);
          i++;
        }
      }
    }
    // run actual forecast
    auto res = forecast_ar(h, newxreg, this->prev_y, this->raw_coef, this->sigma2,
                           this->p, this->intercept, this->drift,
                           this->y.size());
    // scale by sigma2 and take square root
    for( size_t i = 0; i < h; i++ ) {
      res.std_err[i] = sqrt(res.std_err[i] * this->sigma2);
    }
    // rescale forecasts
    if( scalers.size() > 0 ) {
      scalers[0].rescale(res.forecast);
      scalers[0].rescale_w_sd(res.std_err);
    }
    return res;
  };
  const std::vector<scalar_t> get_coef() const { return this->coef; }
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const scalar_t get_sigma2() const { return this->sigma2; }
  const scalar_t get_aic() const { return this->aic; }
  const scalar_t get_bic() const { return this->bic; }
  const scalar_t get_aicc() const { return this->aicc; }
  const bool is_fitted() const { return this->fitted; }
};

#endif
