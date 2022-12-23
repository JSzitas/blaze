#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "initializers.h"
#include "utils/utils.h"
#include "structural_model.h"
#include "xreg.h"
#include "arima_utils.h"
#include "arima_solvers.h"

template <typename U=double> class Arima {
public:
  Arima<U>(){};
  Arima<U>( std::vector<U> & y,
            arima_kind kind,
            std::vector<std::vector<U>> xreg = {{}},
            bool intercept = true,
            bool transform_parameters = true,
            SSinit ss_init = Gardner,
            fitting_method method = ML,
            U kappa = 1000000 ){
    this->y = y;
    // initialize xreg coef and data
    this->xreg = xreg;
    this->intercept = ((kind.d() + kind.D()) == 0 ) && intercept;
    this->transform_parameters = transform_parameters;
    this->ss_init = ss_init;
    this->kind = kind;
    this->method = method;
  };
  void fit(){
    // this should just proceed with fitting, not do things which can be done in
    // the constructor
    // create deltas
    std::vector<U> deltas = make_delta( this->kind.d(), this->kind.period(), this->kind.D());
    // get number of available observations
    int available_n = this->y.size();
    // find na across y
    std::vector<int> na_cases = find_na(y);
    // fit xreg
    if( this->xreg.size() > 0 || this->intercept) {
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
      if( y_d.size() <= xreg_d.size() ) {
        reg_coef = xreg_coef(this->y, this->xreg, this->intercept);
      }
      else {
        reg_coef = xreg_coef(y_d, xreg_d);
      }
      this->reg_coef = reg_coef;
      // find na cases across xreg
      for( int i= 0; i < xreg.size(); i++  ) {
        na_cases = intersect( na_cases, find_na(xreg[i]));
      }
    }
    int missing_cases = na_cases.size();
    available_n -= (deltas.size() + missing_cases);
    // override method to ML if any cases are missing
    if(this->method == CSSML) {
      if(missing_cases > 0) {
        this->method = ML;
      }
    }
    // we have to include xreg in full parameter vector when optimizing -
    // as it can have an impact on the result, it has to be jointly optimized

    // ncond is the number of parameters we are effectively estimating thanks to
    // seasonal parameters
    int ncond = 0;
    if( method == CSS || method == CSSML ) {
      ncond += this->kind.d() + (this->kind.D() * this->kind.period());
      ncond += this->kind.p() + (this->kind.P() * this->kind.period());
    }
    if( ncond <= 0 ) {
      // too few non missing observations - alternatively we can throw
      return;
    }
    // allocate coef vector
    this->coef = std::vector<U>(kind.p() + kind.q() + kind.P() + kind.Q() + reg_coef.size());
    if( method == CSS ) {
      // is using conditional sum of squares, just directly optimize and
      // use hessian 'as-is'
      // print_vector(this->reg_coef.coef);
      const bool is_seasonal = kind.P() + kind.Q();
      if( this->reg_coef.size() > 0 ) {
        if( is_seasonal ) {
          arima_solver_css<true, true>( this->y, this->model, this->reg_coef,
                                        this->xreg, this->kind,
                                        this->coef, ncond );
        }
        else {
          arima_solver_css<true, false>( this->y, this->model, this->reg_coef,
                                         this->xreg, this->kind,
                                         this->coef, ncond );
        }
      }
      else {
        if( is_seasonal ) {
          arima_solver_css<false, true>( this->y, this->model, this->reg_coef,
                                         this->xreg, this->kind,
                                         this->coef, ncond );
        }
        else {
          arima_solver_css<false, false>( this->y, this->model, this->reg_coef,
                                          this->xreg, this->kind,
                                          this->coef, ncond );
        }
      }



      //         res <- optim(init[mask], armaCSS, method = "BFGS",
      //                      hessian = TRUE, control = optim.control)
      //         coef[mask] <- res$par
      //         trarma <- .Call(stats:::C_ARIMA_transPars, coef, arma, FALSE)
      //         mod <- stats:::makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
      //                                  SSinit)
      //         if (ncxreg > 0) {
      //           x <- x - xreg %*% coef[narma + (1L:ncxreg)]
      //         }
      //         val <- .Call(stats:::C_ARIMA_CSS, x, arma, trarma[[1L]], trarma[[2L]],
      //                      as.integer(ncond), TRUE)
      //           sigma2 <- val[[1L]]
      //         var <- solve(res$hessian * n.used)
    }
    else {
      //         if (method == "CSS-ML") {
      //           res <- optim(init[mask], armaCSS, method = "BFGS",
      //                        hessian = FALSE, control = optim.control)
      //           if (res$convergence == 0)
      //             init[mask] <- res$par
      //             if (arma[1L] > 0)
      //               if (!arCheck(init[1L:arma[1L]]))
      //                 stop("non-stationary AR part from CSS")
      //                 if (arma[3L] > 0)
      //                   if (!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
      //                     stop("non-stationary seasonal AR part from CSS")
      //                     ncond <- 0L
      //         }
      //         if (transform.pars) {
      //           init <- .Call(stats:::C_ARIMA_Invtrans, init, arma)
      //           if (arma[2L] > 0) {
      //             ind <- arma[1L] + 1L:arma[2L]
      //             init[ind] <- maInvert(init[ind])
      //           }
      //           if (arma[4L] > 0) {
      //             ind <- sum(arma[1L:3L]) + 1L:arma[4L]
      //             init[ind] <- maInvert(init[ind])
      //           }
      //         }
      //         trarma <- .Call(stats:::C_ARIMA_transPars, init, arma, transform.pars)
      //           mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
      //                            SSinit)
      //           res <- optim(init[mask], armafn, method = "BFGS",
      //                        hessian = TRUE, control = optim.control, trans = as.logical(transform.pars))
      //           if (res$convergence > 0)
      //             warning(gettextf("possible convergence problem: optim gave code = %d",
      //                              res$convergence), domain = NA)
      //             coef[mask] <- res$par
      //             if (transform.pars) {
      //               if (arma[2L] > 0L) {
      //                 ind <- arma[1L] + 1L:arma[2L]
      //                 if (all(mask[ind]))
      //                   coef[ind] <- maInvert(coef[ind])
      //               }
      //               if (arma[4L] > 0L) {
      //                 ind <- sum(arma[1L:3L]) + 1L:arma[4L]
      //                 if (all(mask[ind]))
      //                   coef[ind] <- maInvert(coef[ind])
      //               }
      //               if (any(coef[mask] != res$par)) {
      //                 oldcode <- res$convergence
      //                 res <- optim(coef[mask], armafn, method = "BFGS",
      //                              hessian = TRUE, control = list(maxit = 0L,
      //                                                             parscale = optim.control$parscale), trans = TRUE)
      //                 res$convergence <- oldcode
      //                 coef[mask] <- res$par
      //               }
      //               A <- .Call(stats:::C_ARIMA_Gradtrans, as.double(coef), arma)
      //                 A <- A[mask, mask]
      //               var <- crossprod(A, solve(res$hessian * n.used, A))
      //                 coef <- .Call(stats:::C_ARIMA_undoPars, coef, arma)
      //             }
      //             else var <- solve(res$hessian * n.used)
      //               trarma <- .Call(stats:::C_ARIMA_transPars, coef, arma, FALSE)
      //               mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
      //                                SSinit)
      //               val <- if (ncxreg > 0L) {
      //                 arimaSS(x - xreg %*% coef[narma + (1L:ncxreg)], mod)
      //               } else arimaSS(x, mod)
      //                 sigma2 <- val[[1L]][1L]/n.used
    }
    //       value <- 2 * n.used * res$value + n.used + n.used * log(2 * pi)
    //         aic <- ifelse(method != "CSS", value + 2 * sum(mask) + 2, NA)
    //         resid <- val[[2L]]
  };
  forecast_result<U> predict(int h = 10, std::vector<std::vector<U>> newxreg = {{}}){
     // validate xreg length
     if( newxreg.size() != this->xreg.size() ) {
       return forecast_result<U>(0) ;
     }
     auto res = kalman_forecast(h, this->model);
     auto xreg_adjusted = std::vector<U>(h);
     if( this->reg_coef.size() > 0 ) {
       // get the result of xreg regression
       xreg_adjusted = predict(this->reg_coef, newxreg);
     }
     for(int i=0; i < res.forecast.size(); i++) {
       res.forecast[i] += xreg_adjusted[i];
     }
     // res.forecast = xreg_adjusted + res.forecast;
     for( int i = 0; i < res.se.size(); i++ ) {
       res.se[i] = res.se[i] * this->sigma2;
     }
    return res;
  };
  const std::vector<U> get_coef() const {
    return this->coef;
  }
private:
  std::vector<U> y;
  structural_model<U> model;
  arima_kind kind;
  std::vector<U> coef;
  std::vector<U> residuals;
  std::vector<std::vector<U>> xreg;
  lm_coef<U> reg_coef;
  bool intercept;
  bool transform_parameters;
  SSinit ss_init;
  fitting_method method;
  U sigma2;
};

#endif
