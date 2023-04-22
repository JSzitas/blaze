#ifndef ARIMA_DISPATCH
#define ARIMA_DISPATCH

#include "arima/structures/structures.h"

#include "arima/solvers/ml_solver.h"
#include "arima/solvers/css_solver.h"
#include "arima/solvers/lasso_solver.h"

// it is easier to do dispatch over potential fitting methods here

template<typename scalar_t> scalar_t fit_arima(
             std::vector<scalar_t> &y,
             structural_model<scalar_t> &model,
             const bool intercept,
             const bool drift,
             std::vector<std::vector<scalar_t>> &xreg,
             const arima_kind &kind,
             std::vector<scalar_t> &coef,
             const scalar_t kappa,
             const SSinit ss_init,
             std::vector<scalar_t> &residuals,
             const fitting_method method,
             const bool is_seasonal,
             const bool has_xreg,
             const bool transform_parameters) {
  scalar_t sigma2 = 0;
  if(method == CSSML || CSS) {
 /* this is an ugly tower of specializations, but as far as I can tell,
  * it is the simplest (if ugliest) way to do it
  * to reassure you if you are reading this, all calls are the same,
  * and the only thing that changes are the template arguments
  */
  if (has_xreg) {
    if (is_seasonal) {
      sigma2 = css_solver<true, true, scalar_t>(
        y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
    } else {
      sigma2 = css_solver<true, false, scalar_t>(
        y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
      }
    } else {
      if (is_seasonal) {
        sigma2 = css_solver<false, true, scalar_t>(
          y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
      } else {
        sigma2 = css_solver<false, false, scalar_t>(
          y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
        }
      }
  }
  if( method == ML || method == CSSML) {
    if(transform_parameters) {
      arima_inverse_transform_parameters<
        decltype(coef), scalar_t>(coef, kind);
    }
    /* again, ugly tower, all calls are the same and differ only in template
    * parameters - this is the (sadly) easiest way to do it :/ */
    if (has_xreg) {
      if (is_seasonal) {
        if(transform_parameters) {
          sigma2 = ml_solver<true, true, true, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        } else{
          sigma2 = ml_solver<true, true, false, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        }
      } else {
        if(transform_parameters) {
          sigma2 = ml_solver<true, false, true, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        } else{
          sigma2 = ml_solver<true, false, false, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        }
      }
    } else {
      if (is_seasonal) {
        if(transform_parameters) {
          sigma2 = ml_solver<false, true, true, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        } else{
          sigma2 = ml_solver<false, true, false, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        }
      } else {
        if(transform_parameters) {
          sigma2 = ml_solver<false, false, true, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        } else{
          sigma2 = ml_solver<false, false, false, scalar_t>(
            y, model, intercept, drift, xreg, kind, coef, kappa, ss_init, residuals);
        }
      }
    }
  }
  if(method == LASSO) {
    if (has_xreg) {
      if (is_seasonal) {
        sigma2 = lasso_solver<true, true, scalar_t>(
          y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
      } else {
        sigma2 = lasso_solver<true, false, scalar_t>(
          y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
      }
    } else {
      if (is_seasonal) {
        sigma2 = lasso_solver<false, true, scalar_t>(
          y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
      } else {
        sigma2 = lasso_solver<false, false, scalar_t>(
          y, kind, model, xreg, intercept, drift, coef, kappa, ss_init, residuals);
      }
    }
  }
  return sigma2;
}

#endif
