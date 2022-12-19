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

class ARIMA_CSS_PROBLEM : public FunctionXd {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_CSS_PROBLEM( std::vector<double> &y,
                     arima_kind &kind,
                     structural_model<double> &model,
                     lm_coef<double> & xreg_pars,
                     std::vector<std::vector<double>> & xreg,
                     int n_cond ) : y(y), kind(kind), n_cond(n_cond),
                     model(model), xreg_pars(xreg_pars) {
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
  }
  double operator()( const Eigen::VectorXd &x) const {
    Eigen::VectorXd new_x = x;
    // print_eigvec(new_x);
    Eigen::VectorXd y_temp(y.size());
    for(int i=0; i < y_temp.size(); i++) {
      y_temp[i] = this->y[i];
    }
    if(xreg.cols() > 0) {
      y_temp = y_temp - this->xreg * new_x.tail(new_x.size() - this->arma_pars);
      // do differencing here if necessary - otherwise do differencing when getting
      // the data first - in the constructor
    }

    arima_transform_parameters(new_x, this->kind, false);
    // print_eigvec(new_x);
    // call arima css function
    double res = arima_css_ssq( y_temp, new_x, this->kind, this->n_cond );
    return 0.5 * log(res);
  }
  std::vector<double> y;
  structural_model<double> model;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xreg;
  lm_coef<double> xreg_pars;
  arima_kind kind;
  int n_cond;
  int arma_pars;
};

void arima_solver_css( std::vector<double> &y,
                       structural_model<double> &model,
                       lm_coef<double> xreg_coef,
                       std::vector<std::vector<double>> xreg,
                       arima_kind kind,
                       int n_cond ) {
  // define solver
  using Solver = cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM>;

  // initialize a solver object
  Solver solver;

  // initialize problem
  ARIMA_CSS_PROBLEM css_arima_problem( y, kind, model, xreg_coef,
                                       xreg, n_cond);
  // initialize function input
  // first determine size >>
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
  // print_eigvec(x);
  // arima_transform_parameters(x, kind, false);
  // print_eigvec(x);

  // auto init_solution = css_arima_problem.Eval(x);
  // print_eigvec(init_solution.x);
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // print solution, solver state and estimated coef
  std::vector<double> result(x.size());
  for( int i=0; i < result.size(); i++) {
    result[i] = solution.x[i];
  }
  print_vector(result);
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
