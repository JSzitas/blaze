#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "utils/utils.h"

#include "arima/structures/fitting_method.h"
#include "arima/structures/structural_model.h"

#include "arima/utils/xreg.h"

// #include "arima/solvers/arima_css_solver.h"
#include "arima/solvers/arima_css_optim_grad.h"
#include "arima/solvers/arima_ml_solver.h"

#include "arima/utils/checks.h"


template <typename U = double, typename Scaler = StandardScaler<U>> class Arima {
private:
  // data
  std::vector<U> y;
  structural_model<U> model;
  arima_kind kind;
  std::vector<std::vector<U>> xreg;
  bool intercept, transform_parameters;
  SSinit ss_init;
  fitting_method method;
  U kappa;
  // estimated during fitting
  std::vector<U> coef, residuals, reg_coef;
  U sigma2;
  std::vector<Scaler> scalers;
  U aic;
  std::array<bool, 2> ar_stationary;
  bool fitted;
public:
  Arima<U, Scaler>(){};
  Arima<U, Scaler>(
      const std::vector<U> &y, const arima_kind kind,
      const std::vector<std::vector<U>> xreg = {{}}, const bool intercept = true,
      const bool transform_parameters = true, const SSinit ss_init = Gardner,
      const fitting_method method = ML, const U kappa = 1000000,
      const bool standardize = true) : y(y), kind(kind), xreg(xreg),
      intercept(((kind.d() + kind.D()) == 0) && intercept),
      transform_parameters(transform_parameters),
      ss_init(ss_init), method(method), kappa(kappa) {
    this->scalers = std::vector<Scaler>( standardize * (1 + xreg.size()) );
    this->model = structural_model<U>();
    this->fitted = false;
    this->ar_stationary = {true,true};
  };
  void fit() {
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
    std::vector<U> reg_coef( this->xreg.size(), this->intercept );
    // we have to include xreg in full parameter vector when optimizing -
    // as it can have an impact on the result, it has to be jointly optimized
    // ncond is the number of parameters we are effectively estimating thanks to
    // seasonal parameters
    // allocate coef vector
    const size_t arma_coef_size = this->kind.p() + this->kind.q() + this->kind.P() + this->kind.Q();
    this->coef = std::vector<U>(arma_coef_size + reg_coef.size(), 0);
    // load any estimated xreg into coef
    // fit xreg
    if (this->xreg.size() > 0 || this->intercept) {
      std::vector<U> y_d;
      std::vector<std::vector<U>> xreg_d;
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
      reg_coef = xreg_coef(this->y, this->xreg, this->intercept);
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
    if (this->method != ML) {
      /* this is an ugly tower of specializations, but as far as I can tell,
       * it is the simplest (if ugliest) way to do it
       * to reassure you if you are reading this, all calls are the same,
       * and the only thing that changes are the template arguments
       */
      if (this->reg_coef.size() > 0) {
        if (is_seasonal) {
          this->sigma2 = arima_solver_css<true, true>(
            this->y, this->kind, this->model,
            this->xreg, this->intercept, this->coef,
            this->kappa, this->ss_init);
        } else {
          arima_solver_css<true, false>(
              this->y, this->kind, this->model,
              this->xreg, this->intercept, this->coef,
              this->kappa, this->ss_init);
        }
      } else {
        if (is_seasonal) {
          arima_solver_css<false, true>(
              this->y, this->kind, this->model,
              this->xreg, this->intercept, this->coef,
              this->kappa, this->ss_init);
        } else {
          arima_solver_css<false, false>(
              this->y, this->kind, this->model,
              this->xreg, this->intercept, this->coef,
              this->kappa, this->ss_init);
        }
      }
    }
    if( this->method == CSSML) {
      //perform checks on AR coefficients following CSS fit
      // diagnostic struct that carries codes for unstable fits?
      // would allow us to have nothrow all over these :)
      this->ar_stationary = check_all_ar(this->coef, this->kind);
    }
    if( this->method == ML || this->method == CSSML) {
      if(this->transform_parameters) {
        arima_inverse_transform_parameters(this->coef, this->kind);
      }
      /* again, ugly tower, all calls are the same and differ only in template
       * parameters - this is the (sadly) easiest way to do it :/
       */
      if (this->reg_coef.size() > 0) {
        if (is_seasonal) {
          if(this->transform_parameters) {
            this->sigma2 = arima_solver_ml<true, true, true>(
              this->y, this->model, this->intercept,
              this->xreg, this->kind, this->coef,
              this->kappa, this->ss_init);
          }
          else{
            this->sigma2 = arima_solver_ml<true, true, false>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          }
        } else {
          if(this->transform_parameters) {
            this->sigma2 = arima_solver_ml<true, false, true>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          }
          else{
            this->sigma2 = arima_solver_ml<true, false, false>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          }
        }
      } else {
        if (is_seasonal) {
          if(this->transform_parameters) {
            this->sigma2 = arima_solver_ml<false, true, true>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          } else{
            this->sigma2 = arima_solver_ml<false, true, false>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          }
        } else {
          if(this->transform_parameters) {
            this->sigma2 = arima_solver_ml<false, false, true>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          } else{
            this->sigma2 = arima_solver_ml<false, false, false>(
                this->y, this->model, this->intercept,
                this->xreg, this->kind, this->coef,
                this->kappa, this->ss_init);
          }
        }
      }
    }
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
      constexpr double one_p_log_twopi = 1.0 + 1.837877;
      this->aic = available_n * (log(sigma_2) + one_p_log_twopi);
    }
    // invert scaling
    if( this->scalers.size() > 0 ) {
      // first scaler used for target
      this->scalers[0].rescale(this->y);
      size_t i = 1;
      for( auto & xreg_val:this->xreg ) {
        U temp = this->scalers[i].rescale_val(coef[arma_coef_size+i-1]);
        this->coef[arma_coef_size+i-1] = temp;
        i++;
      }
      if( this->intercept ) {
        // the intercept also has to be rescaled to have any meaning
        U temp = scalers[0].rescale_val(this->coef.back());
        this->coef.back() = temp;
      }
    }
    this->fitted = true;
  };
  forecast_result<U> forecast(const size_t h = 10,
                              std::vector<std::vector<U>> newxreg = {{}}) {
    // validate xreg length
    if (!this->fitted || newxreg.size() != this->xreg.size())
      return forecast_result<U>(0);
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
      auto xreg_adjusted = predict(h, this->reg_coef, this->intercept, newxreg);
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
  const structural_model<U> get_structural_model() const { return this->model; }
  const std::vector<U> get_coef() const { return this->coef; }
  const bool is_fitted() const { return this->fitted; }
};

#endif
