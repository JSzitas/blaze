#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "initializers.h"
#include "utils/utils.h"
#include "structural_model.h"
#include "xreg.h"

// defines the arima structure
struct arima_kind{
   arima_kind( int p, int d, int q, int P, int D, int Q, int s_period ) :
    p(p), diff_d(d), q(q), P(P), seas_diff_D(D), Q(Q), s_period(s_period){}
  int d(){
    return this->diff_d;
  }
  int D(){
    return this->seas_diff_D;
  }
  int period() {
    return this->s_period;
  }
private:
  int p, diff_d, q, P, seas_diff_D, Q, s_period;
};

template <typename U=double> class Arima {
  Arima<U>(std::vector<U> & y,
           arima_kind kind,
           std::vector<std::vector<U>> xreg = {{}},
           bool intercept = true,
           bool transform_parameters = true,
           SSinit ss_init = Gardner,
           U kappa = 1000000 ){
    this->y = y;
    this->xreg = xreg;
    this->intercept = intercept;
    this->transform_parameters = transform_parameters;
    this->ss_init = ss_init;
    this->kind = kind;
    // this->arma_structure = std::vector<int>{ order[0], order[2], seas_order[0],
    //                                          seas_order[2], seasonal_period,
    //                                          order[1], seas_order[1]};
    // for(int i=0; i < 4; i++) {
    //   this->n_arma_coef += this->arma_structure[i];
    // }

  };
  void fit(){
    // this should just proceed with fitting, not do things which can be done in
    // the constructor

    // take a copy so we do not modify the underlying data - maybe we
    // change this later
    std::vector<U> y_fit = this->y;
    // fit xreg
    if( this->xreg.size() > 0 ) {
      std::vector<U> y_d;
      std::vector<std::vector<U>> xreg_d;
      // if we have any differences
      if( this->kind.d() > 0 ) {
        y_d = diff(this->y, 1, this->kind.d());
        xreg_d = diff(this->xreg, 1, this->kind.d());
      }
      // seasonal differences
      if( this->kind.period() > 1 && this->kind.D() > 0  ) {
        y_d = diff(this->y, this->kind.period(), this->kind.D());
        xreg_d = diff(this->xreg, this->kind.period(), this->kind.D());
      }
      lm_coef<U> reg_coef(xreg_d.size(), this->intercept);
      // fit coefficients and adjust y for fitted coefficients -
      // the original R code does this repeatedly, but it can be done only once
      // - the fitted effects from external regressors are never refitted
      if( this->y_d <= xreg_d.size() ) {
        reg_coef = xreg_coef(this->y, this->xreg, this->intercept);
        y_fit -= predict(reg_coef, xreg);
      }
      else {
        reg_coef = xreg_coef(y_d, xreg_d);
        y_fit -= predict(reg_coef, xreg);
      }
      this->reg_coef = reg_coef;

      //         isna <- is.na(x) | apply(xreg, 1L, anyNA)
      //           n.used <- sum(!isna) - length(Delta)
      //           init0 <- c(init0, coef(fit))
      //           ses <- summary(fit)$coefficients[, 2L]
      //         parscale <- c(parscale, 10 * ses)
      // }
    }



  };
  forecast_result<U> predict(int h = 10, std::vector<std::vector<U>> newxreg = {{}}){
     // validate xreg length
     if( newxreg.size() != this->xreg.size() ) {
       return forecast_result<U>(0) ;
     }
     auto res = kalman_forecast(h, this-> structural_arma_model);
     auto xreg_adjusted = std::vector<U>(h);
     if( this->reg_coef.size() > 0 ) {
       // get the result of xreg regression
       xreg_adjusted = predict(this->reg_coef, newxreg);
     }
     res.forecast = xreg_adjusted + res.forecast;
     for( int i = 0; i < res.se.size(); i++ ) {
       res.se[i] = res.se[i] * sigma2;
     }
    return res;
  };
  std::vector<U> y;
  arima_kind kind;
  std::vector<U> residuals;
  std::vector<std::vector<U>> xreg;
  lm_coef<U> reg_coef;
  bool intercept;
  bool transform_parameters;
  SSinit ss_init;
  U sigma2;
};

#endif
