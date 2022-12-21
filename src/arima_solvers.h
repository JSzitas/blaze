#ifndef ARIMA_SOLVER
#define ARIMA_SOLVER

#include "utils/utils.h"
// included mainly for isnan()
#include <math.h>

#include "arima_css_likelihood.h"
#include "arima_utils.h"
#include "structural_model.h"
#include "xreg.h"
// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

enum fitting_method{
  CSS = 1,
  CSSML = 2,
  ML = 3
};

using FunctionXd = cppoptlib::function::Function<double>;

template <const bool has_xreg, const bool seasonal>class ARIMA_CSS_PROBLEM : public FunctionXd {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_CSS_PROBLEM( std::vector<double> &y,
                     arima_kind &kind,
                     structural_model<double> &model,
                     lm_coef<double> & xreg_pars,
                     std::vector<std::vector<double>> & xreg,
                     int n_cond ) : kind(kind), n_cond(n_cond),
                     xreg_pars(xreg_pars) {
    // initialize an xreg matrix
    int n = y.size();
    std::vector<double> _xreg = flatten_vec(xreg);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> new_mat = Eigen::Map<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    >(_xreg.data(), n, xreg.size()/n);
    if( xreg_pars.has_intercept() ) {
      new_mat.conservativeResize(Eigen::NoChange, new_mat.cols()+1);
      new_mat.col(new_mat.cols()-1) = Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(n, 1, 1);
    }
    this->xreg = new_mat;
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    Eigen::VectorXd y_temp(n);
    for( int i=0; i < n; i++ ) {
      y_temp[i] = y[i];
    }
    if constexpr(!has_xreg) {
      // if you have no xreg, you can pull differencing out and only do it once
      // rather than do it repeatedly at every step
      // regular differencing
      for (int i = 0; i < kind.d(); i++) {
        for (int l = n - 1; l > 0; l--) {
          y_temp[l] -= y_temp[l - 1];
        }
      }
      // seasonal differencing
      int ns = kind.period();
      for (int i = 0; i < kind.D(); i++) {
        for (int l = n - 1; l >= ns; l--) {
          y_temp[l] -= y_temp[l - ns];
        }
      }
    }
    this->y = y;
    this->y_temp = y_temp;
    this->n = n;
    this->xreg = new_mat;
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    // preallocate new_x
    this->new_x = Eigen::VectorXd(
      kind.p() + (kind.P() * kind.period()) +
      kind.q() + (kind.Q() *kind.period()) +
      this->xreg.cols());
    // preallocate model residuals
    this->residual = std::vector<double>(n);
  }
  double operator()( const Eigen::VectorXd &x) {
    for( int i=0; i < x.size(); i++ ) {
      this->new_x[i] = x[i];
    }

    if constexpr(has_xreg) {
      // refresh y_temp and load it with original y data
      for(int i=0; i < this->n; i++) {
        this->y_temp[i] = this->y[i];
      }
      if(xreg.cols() > 0) {
        this->y_temp = this->y_temp - this->xreg * x.tail(x.size() - this->arma_pars);
        // do differencing here
        for (int i = 0; i < kind.d(); i++) {
          for (int l = n - 1; l > 0; l--) {
            this->y_temp[l] -= this->y_temp[l - 1];
          }
        }
        // seasonal differencing
        int ns = kind.period();
        for (int i = 0; i < kind.D(); i++) {
          for (int l = n - 1; l >= ns; l--) {
            this->y_temp[l] -= this->y_temp[l - ns];
          }
        }
      }
    }
    arima_transform_parameters<seasonal, false>(this->new_x, this->kind);
    // call arima css function
    double res = arima_css_ssq( this->y_temp, this->new_x, this->kind, this->n_cond, this->residual );
    return 0.5 * log(res);
  }
  std::vector<double> y;
  Eigen::VectorXd y_temp;
  Eigen::VectorXd new_x;
  // structural_model<double> model;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xreg;
  lm_coef<double> xreg_pars;
  std::vector<double> residual;
  arima_kind kind;
  int n_cond;
  int arma_pars;
  int n;
};

template <const bool has_xreg,
          const bool seasonal> void arima_solver_css(
              std::vector<double> &y,
              structural_model<double> &model,
              lm_coef<double> xreg_coef,
              std::vector<std::vector<double>> xreg,
              arima_kind kind,
              int n_cond ) {

  auto vec_size = kind.p() + kind.q() + kind.P() + kind.Q() + xreg_coef.size();
  auto arma_size = kind.p() + kind.q() + kind.P() + kind.Q();
  Eigen::VectorXd x(vec_size);
  for( auto &val:x ) {
    val = 0;
  }
  // initialize to all zeroes except for xreg
  for(int i = arma_size; i < vec_size; i++) {
    x[i] = xreg_coef.coef[i-arma_size];
  }
  std::vector<double> result(x.size());
  using Solver = cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal>>;
  Solver solver;
  ARIMA_CSS_PROBLEM<has_xreg, seasonal> css_arima_problem( y, kind, model, xreg_coef,
                                                           xreg, n_cond);
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  for( int i=0; i < result.size(); i++) {
    result[i] = solution.x[i];
  }
  // print_vector(result);
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








#endif
