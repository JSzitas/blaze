#ifndef ARIMA_WRAPPER
#define ARIMA_WRAPPER

#include "arima/solvers/arima_css_solver.h"
#include "arima/solvers/initializers.h"
#include "arima/structures/structural_model.h"
#include "arima/utils/arima_utils.h"
#include "arima/utils/xreg.h"
#include "utils/utils.h"

template <typename U = double> class Arima {
public:
  Arima<U>(){};
  Arima<U>(std::vector<U> &y, arima_kind kind,
           std::vector<std::vector<U>> xreg = {{}}, bool intercept = true,
           bool transform_parameters = true, SSinit ss_init = Gardner,
           fitting_method method = ML, U kappa = 1000000) {
    this->y = y;
    // initialize xreg coef and data
    this->xreg = xreg;
    this->intercept = ((kind.d() + kind.D()) == 0) && intercept;
    this->transform_parameters = transform_parameters;
    this->ss_init = ss_init;
    this->kind = kind;
    this->method = method;
    this->kappa = kappa;
    this->model = structural_model<U>();
    this->fitted = false;
  };
  void fit() {
    // this should just proceed with fitting, not do things which can be done in
    // the constructor
    // create deltas
    std::vector<U> deltas =
        make_delta(this->kind.d(), this->kind.period(), this->kind.D());
    // get number of available observations
    size_t available_n = this->y.size();
    // find na across y
    std::vector<size_t> na_cases = find_na(this->y);
    // initialize xreg
    lm_coef<U> reg_coef( this->xreg.size(), this->intercept );;
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
      // fit coefficients and adjust y for fitted coefficients -
      // the original R code does this repeatedly, but it can be done only once
      // - the fitted effects from external regressors are never refitted
      if (y_d.size() <= xreg_d.size()) {
        reg_coef = xreg_coef(this->y, this->xreg, this->intercept);
      } else {
        reg_coef = xreg_coef(y_d, xreg_d);
      }
      // find na cases across xreg
      for (size_t i = 0; i < this->xreg.size(); i++) {
        na_cases = intersect(na_cases, find_na(this->xreg[i]));
      }
    }
    // store regression coefficients (if any) in this object
    this->reg_coef = reg_coef;
    // adjust CSS for missing cases
    size_t missing_cases = na_cases.size();
    available_n -= (deltas.size() + missing_cases);
    // override method to ML if any cases are missing
    if (this->method == CSSML) {
      if (missing_cases > 0) {
        this->method = ML;
      }
    }
    // we have to include xreg in full parameter vector when optimizing -
    // as it can have an impact on the result, it has to be jointly optimized
    // ncond is the number of parameters we are effectively estimating thanks to
    // seasonal parameters
    int ncond = 0;
    if (this->method == CSS || this->method == CSSML) {
      ncond += this->kind.d() + (this->kind.D() * this->kind.period());
      ncond += this->kind.p() + (this->kind.P() * this->kind.period());
    }
    if (ncond <= 0) {
      // too few non missing observations - alternatively we can throw
      return;
    }
    // allocate coef vector
    const int arma_coef_size = this->kind.p() + this->kind.q() + this->kind.P() + this->kind.Q();
    this->coef = std::vector<U>(arma_coef_size + reg_coef.size());
    if (this->method == CSS) {
      // is using conditional sum of squares, just directly optimize and
      const bool is_seasonal = this->kind.P() + this->kind.Q();
      const bool has_xreg = this->reg_coef.size() > 0;
      if (this->reg_coef.size() > 0) {
        if (is_seasonal) {
          arima_solver_css<true, true>(this->y, this->model, this->reg_coef,
                                       this->xreg, this->kind, this->coef,
                                       deltas, ncond, available_n, this->kappa,
                                       this->ss_init, this->sigma2);
        } else {
          arima_solver_css<true, false>(this->y, this->model, this->reg_coef,
                                        this->xreg, this->kind, this->coef,
                                        deltas, ncond, available_n, this->kappa,
                                        this->ss_init, this->sigma2);
        }
      } else {
        if (is_seasonal) {
          arima_solver_css<false, true>(this->y, this->model, this->reg_coef,
                                        this->xreg, this->kind, this->coef,
                                        deltas, ncond, available_n, this->kappa,
                                        this->ss_init, this->sigma2);
        } else {
          arima_solver_css<false, false>(
              this->y, this->model, this->reg_coef, this->xreg, this->kind,
              this->coef, deltas, ncond, available_n, this->kappa,
              this->ss_init, this->sigma2);
        }
      }
      // load xreg coefficients from coef as necessary
      for (size_t i = arma_coef_size; i < this->coef.size(); i++) {
        this->reg_coef[i - arma_coef_size] = this->coef[i];
      }
    }
    if( this->method == CSSML) {
      //perform checks on AR coefficients
      // if( !ar_check() )
      //             if (!arCheck(init[1L:arma[1L]]))
      //               stop("non-stationary AR part from CSS")
      //               if (arma[3L] > 0)
      //                 if (!arCheck(init[sum(arma[1L:2L]) + 1L:arma[3L]]))
      //                   stop("non-stationary seasonal AR part from CSS")
      //                   ncond <- 0L
    }
    if( this->method == ML || this->method == CSSML) {
            // if (this->transform_parameters) {
            //   init <- .Call(stats:::C_ARIMA_Invtrans, init, arma)
            //   if (this->kind.q()) {
            //     ind <- arma[1L] + 1L:arma[2L]
            //     init[ind] <- maInvert(init[ind])
            //   }
            //   if (this->kind.Q()) {
            //     ind <- sum(arma[1L:3L]) + 1L:arma[4L]
            //     init[ind] <- maInvert(init[ind])
            //   }
            // }
            // trarma <- .Call(stats:::C_ARIMA_transPars, init, arma,
            // transform.pars)
            //   mod <- makeARIMA(trarma[[1L]], trarma[[2L]], Delta, kappa,
            //                    SSinit)
            //   res <- optim(init[mask], armafn, method = "BFGS",
            //                hessian = TRUE, control = optim.control, trans =
            //                as.logical(transform.pars))
            //   if (res$convergence > 0)
            //     warning(gettextf("possible convergence problem: optim gave
            //     code = %d",
            //                      res$convergence), domain = NA)
            //     coef[mask] <- res$par
            //     if (transform.pars) {
            //       if (arma[2L] > 0L) {
            //         ind <- arma[1L] + 1L:arma[2L]
            //         if (all(mask[ind]))
            //           coef[ind] <- maInvert(coef[ind])
            //       }
            //       if (arma[4L] > 0L) {
            //         ind <- sum(arma[1L:3L]) + 1L:arma[4L]
            //         if (all(mask[ind]))
            //           coef[ind] <- maInvert(coef[ind])
            //       }
            //       if (any(coef[mask] != res$par)) {
            //         oldcode <- res$convergence
            //         res <- optim(coef[mask], armafn, method = "BFGS",
            //                      hessian = TRUE, control = list(maxit = 0L,
            //                                                     parscale =
            //                                                     optim.control$parscale),
            //                                                     trans =
            //                                                     TRUE)
            //         res$convergence <- oldcode
            //         coef[mask] <- res$par
            //       }
            //       A <- .Call(stats:::C_ARIMA_Gradtrans, as.double(coef),
            //       arma)
            //         A <- A[mask, mask]
            //       var <- crossprod(A, solve(res$hessian * n.used, A))
            //         coef <- .Call(stats:::C_ARIMA_undoPars, coef, arma)
            //     }
            //     else var <- solve(res$hessian * n.used)
            //       trarma <- .Call(stats:::C_ARIMA_transPars, coef, arma,
            //       FALSE) mod <- makeARIMA(trarma[[1L]], trarma[[2L]],
            //       Delta, kappa,
            //                        SSinit)
            //       val <- if (ncxreg > 0L) {
            //         arimaSS(x - xreg %*% coef[narma + (1L:ncxreg)], mod)
            //       } else arimaSS(x, mod)
            //         sigma2 <- val[[1L]][1L]/n.used
    }
    if (this->method != CSS) {
      this->aic = std::nan("");
    } else {
      // 1.837877 is equal to log(2*pi) - log is not standard compliant in a
      // constexpr so we must expand the expression manually, sadly
      constexpr double one_p_log_twopi = 1.0 + 1.837877;
      this->aic = available_n * (log(this->sigma2) + one_p_log_twopi);
    }
    this->fitted = true;
  };
  forecast_result<U> forecast(const size_t h = 10,
                              std::vector<std::vector<U>> newxreg = {{}}) {
    // validate xreg length
    if (!this->fitted || newxreg.size() != this->xreg.size()) {
      return forecast_result<U>(0);
    }
    // rescale sigma2 - this is necessary because the original sigma2 is
    // a conditional sum of squares
    double rescaled_sigma2 = exp(0.5 * log(this->sigma2));
    auto res = kalman_forecast(h, this->model, rescaled_sigma2);
    if( reg_coef.has_intercept() ) {
      // handle intercept
      for (size_t i = 0; i < h; i++) {
        res.forecast[i] += this->reg_coef.get_intercept();
      }
    }
    if(newxreg.size() > 0) {
      auto xreg_adjusted = predict(h, this->reg_coef, newxreg);
      for (size_t i = 0; i < h; i++) {
        res.forecast[i] += xreg_adjusted[i];
      }
    }
    return res;
  };
  const structural_model<U> get_structural_model() const { return this->model; }
  const std::vector<U> get_coef() const { return this->coef; }
  const bool is_fitted() const { return this->fitted; }

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
  U kappa;
  U aic;
  bool fitted;
};

#endif
