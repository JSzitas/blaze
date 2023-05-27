#ifndef BLAZE_ARIMA
#define BLAZE_ARIMA

#include "utils/utils.h"
#include "utils/xreg.h"

#include "arima/structures/structures.h"
#include "arima/solvers/dispatch.h"

template <typename scalar_t, typename Scaler = StandardScaler<scalar_t>> class Arima {
private:
  // data
  std::vector<scalar_t> y;
  structural_model<scalar_t> model;
  arima_kind kind;
  std::vector<std::vector<scalar_t>> xreg;
  bool intercept, drift, transform_parameters;
  SSinit ss_init;
  fitting_method method;
  scalar_t kappa;
  // estimated during fitting
  std::vector<scalar_t> coef, residuals, fitted_vals, reg_coef;
  scalar_t sigma2;
  std::vector<Scaler> scalers;
  scalar_t aic;
  std::array<bool, 2> ar_stationary;
  bool fitted;
public:
  Arima<scalar_t, Scaler>(){};
  Arima<scalar_t, Scaler>(
      const std::vector<scalar_t> &y,
      const arima_kind &kind,
      const std::vector<std::vector<scalar_t>> &xreg = {{}},
      const bool intercept = true, const bool drift = true,
      const bool transform_parameters = true,
      const SSinit ss_init = SSinit::Gardner,
      const fitting_method method = ML, const scalar_t kappa = 1000000,
      const bool standardize = true) : y(y), kind(kind), xreg(xreg),
      intercept(((kind.d() + kind.D()) == 0) && intercept),
      drift(((kind.d() + kind.D()) == 1) && drift),
      transform_parameters(transform_parameters),
      ss_init(ss_init), method(method), kappa(kappa) {
    this->scalers = std::vector<Scaler>( standardize * (1 + xreg.size()) );
    this->model = structural_model<scalar_t>();
    this->fitted = false;
    this->ar_stationary = {true,true};
    this->residuals = std::vector<scalar_t>(y.size(), 0.0);
    this->fitted_vals = std::vector<scalar_t>(y.size(), 0.0);
  };
  void fit() {
    // if we are asked to fit a model with p,q, P,Q == 0, we simply throw 
    if( (kind.p() + kind.q() + kind.P() + kind.Q()) == 0) {
      std::cout << "You specified a model with no parameters (p,q,P,Q all equal" <<
                   " to zero) - we will not fit a model in this case. Please" <<
                   " use the .respecify() method and provide a valid model." <<
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
    // get number of available observations
    size_t available_n = this->y.size();
    // find na across y
    std::vector<size_t> na_cases = find_na(this->y);
    // initialize xreg
    std::vector<scalar_t> reg_coef( this->xreg.size() + this->intercept + this->drift );
    // we have to include xreg in full parameter vector when optimizing -
    // as it can have an impact on the result, it has to be jointly optimized
    // ncond is the number of parameters we are effectively estimating thanks to
    // seasonal parameters
    // allocate coef vector
    const size_t arma_coef_size = this->kind.p() + this->kind.q() + this->kind.P() + this->kind.Q();
    this->coef = std::vector<scalar_t>(arma_coef_size + reg_coef.size(), 0);
    // load any estimated xreg into coef
    // fit xreg
    if (this->xreg.size() > 0 || this->intercept) {
      std::vector<scalar_t> y_d;
      std::vector<std::vector<scalar_t>> xreg_d;
      // if we have any differences
      if (this->kind.d() > 0) {
        y_d = diff(this->y, 1, this->kind.d());
        xreg_d = diff(this->xreg, 1, this->kind.d());
      }
      // seasonal differences
      if (this->kind.period() > 1 && this->kind.D() > 0) {
        y_d = diff(this->y, this->kind.period(), this->kind.D());
        xreg_d = diff(this->xreg, this->kind.period(), this->kind.D());
      }
      // fit coefficients to initialize fitting
      if( kind.d() || kind.D() ) {
        reg_coef = xreg_coef(y_d, xreg_d, this->intercept, this->drift);
      }
      else {
        reg_coef = xreg_coef(this->y, this->xreg, this->intercept, this->drift);
      }
      // find na cases across xreg
      for (size_t i = 0; i < this->xreg.size(); i++) {
        na_cases = intersect(na_cases, find_na(this->xreg[i]));
      }
      // fill in coef using estimated regression coefficients
      for (size_t i = arma_coef_size; i < coef.size(); i++)
        coef[i] = reg_coef[i - arma_coef_size];
    }
    // store regression coefficients (if any) in this object
    this->reg_coef = reg_coef;
    // adjust CSS for missing cases
    size_t missing_cases = na_cases.size();
    available_n -= ((this->kind.d() + 1) +
      (this->kind.period() * this->kind.D()) +
      missing_cases);
    // override method to ML if any cases are missing
    if (this->method == CSSML) {
      if (missing_cases > 0) {
        this->method = ML;
      }
    }
    // yet unused flag for optimizer - maybe we use this in the future
    // bool optimization_failed = false;
    const bool is_seasonal = this->kind.P() + this->kind.Q();
    const bool has_xreg = this->reg_coef.size() > 0;
    this->sigma2 = fit_arima( this->y, this->model, this->intercept,
                              this->drift, this->xreg, this->kind, this->coef,
                              this->kappa, this->ss_init, this->residuals,
                              this->method, is_seasonal, this->xreg.size() > 0,
                              this->transform_parameters);
    // load xreg coefficients from coef as necessary
    for (size_t i = arma_coef_size; i < this->coef.size(); i++) {
      this->reg_coef[i - arma_coef_size] = this->coef[i];
    }
    if (this->method == CSS) {
      this->aic = std::nan("");
    } else {
      // we have to rescale the sigma to be on the same scale as original data
      auto sigma_2 = this->sigma2;
      if( scalers.size() > 0 ) {
        sigma_2 = scalers[0].rescale_val_w_mean(this->sigma2);
      }
      // 1.837877 is equal to log(2*pi) - log is not standard compliant in a
      // constexpr so we must expand the expression manually, sadly
      constexpr scalar_t one_p_log_twopi = 1.0 + 1.837877;
      this->aic = available_n * (log(sigma_2) + one_p_log_twopi);
    }
    // invert scaling - this is probably redundant, we do not care for estimated
    // coefficients too much (as the main goal of the package is forecasting)
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0].rescale(this->y);
      for(size_t i = 1; i < this->scalers.size(); i++) {
        scalar_t temp = this->scalers[i].rescale_val(coef[arma_coef_size+i-1]);
        this->coef[arma_coef_size+i-1] = temp;
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
    for(size_t i = 0; i < this->y.size(); i++) {
      this->fitted_vals[i] = this->y[i] - this->residuals[i];
    }
  };
  void respecify(const arima_kind &kind) {
    this->kind = kind;
  }
  forecast_result<scalar_t> forecast(
      const size_t h = 10,
      std::vector<std::vector<scalar_t>> newxreg = {{}}) {
    // validate xreg length
    if (!this->fitted || newxreg.size() != this->xreg.size())
      return forecast_result<scalar_t>(0);
    // otherwise run forecast
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
                                   this->drift, newxreg, this->y.size());
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
  const structural_model<scalar_t> get_structural_model() const { return this->model; }
  const std::vector<scalar_t> get_coef() const { return this->coef; }
  const std::vector<scalar_t> get_residuals() const { return this->residuals; }
  const std::vector<scalar_t> get_fitted() const { return this->fitted_vals; }
  const bool is_fitted() const { return this->fitted; }
};

#endif
