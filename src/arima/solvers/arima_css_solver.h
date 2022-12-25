#ifndef ARIMA_SOLVER
#define ARIMA_SOLVER

#include "utils/utils.h"

#include "arima/structures/structural_model.h"

#include "arima/solvers/arima_css_likelihood.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/arima_utils.h"
#include "arima/utils/xreg.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

enum fitting_method { CSS = 1, CSSML = 2, ML = 3 };

using FunctionXd = cppoptlib::function::Function<double>;

template <const bool has_xreg, const bool seasonal>
class ARIMA_CSS_PROBLEM : public FunctionXd {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_CSS_PROBLEM(std::vector<double> &y, const arima_kind &kind,
                    lm_coef<double> &xreg_pars,
                    std::vector<std::vector<double>> &xreg, int n_cond)
      : kind(kind), n_cond(n_cond), xreg_pars(xreg_pars) {
    // initialize an xreg matrix
    int n = y.size();
    std::vector<double> _xreg = flatten_vec(xreg);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> new_mat =
        Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(
            _xreg.data(), n, xreg.size() / n);
    if (xreg_pars.has_intercept()) {
      new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
      new_mat.col(new_mat.cols() - 1) =
          Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(n, 1, 1);
    }
    this->xreg = new_mat;
    this->arma_pars = kind.p() + kind.P() + kind.q() + kind.Q();
    Eigen::VectorXd y_temp(n);
    for (int i = 0; i < n; i++) {
      y_temp[i] = y[i];
    }
    if constexpr (!has_xreg) {
      // if you have no intercept, you can do differencing
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
    // pre-allocate new_x
    this->new_x =
        Eigen::VectorXd(kind.p() + (kind.P() * kind.period()) + kind.q() +
                        (kind.Q() * kind.period()) + this->xreg.cols());
    // pre-allocate model residuals
    this->residual = std::vector<double>(n);
    // pre-allocate transformation helper vector
    this->transform_temp_phi =
        std::vector<double>(kind.p() + (kind.P() * kind.period()));
    this->transform_temp_theta =
        std::vector<double>(kind.q() + (kind.Q() * kind.period()));
  }
  double operator()(const Eigen::VectorXd &x) {
    for (int i = 0; i < x.size(); i++) {
      this->new_x[i] = x[i];
    }
    if constexpr (has_xreg) {
      // refresh y_temp and load it with original y data
      for (int i = 0; i < this->n; i++) {
        this->y_temp[i] = this->y[i];
      }
      this->y_temp =
          this->y_temp - this->xreg * x.tail(x.size() - this->arma_pars);
      // do differencing here
      if (!xreg_pars.has_intercept()) {
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
    arima_transform_parameters<seasonal, false>(this->new_x, this->kind,
                                                this->transform_temp_phi,
                                                this->transform_temp_theta);
    // call arima css function
    double res = arima_css_ssq(this->y_temp, this->new_x, this->kind,
                               this->n_cond, this->residual);
    return 0.5 * log(res);
  }
  std::vector<double> y;
  Eigen::VectorXd y_temp;
  Eigen::VectorXd new_x;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> xreg;
  lm_coef<double> xreg_pars;
  std::vector<double> residual;
  arima_kind kind;
  std::vector<double> transform_temp_phi;
  std::vector<double> transform_temp_theta;
  int n_cond;
  int arma_pars;
  int n;
};

template <const bool has_xreg, const bool seasonal>
void arima_solver_css(std::vector<double> &y, structural_model<double> &model,
                      lm_coef<double> xreg_coef,
                      std::vector<std::vector<double>> xreg,
                      const arima_kind &kind, std::vector<double> &coef,
                      std::vector<double> &delta, const int n_cond,
                      const int n_available, const double kappa,
                      const SSinit ss_init, double &sigma2) {

  auto vec_size = kind.p() + kind.q() + kind.P() + kind.Q() + xreg_coef.size();
  auto arma_size = kind.p() + kind.q() + kind.P() + kind.Q();
  Eigen::VectorXd x(vec_size);
  for (auto &val : x) {
    val = 0;
  }
  // initialize to all zeroes except for xreg
  for (int i = arma_size; i < vec_size; i++) {
    x[i] = xreg_coef.coef[i - arma_size];
  }
  // initialize solver
  using Solver = cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal>>;
  Solver solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal> css_arima_problem(y, kind, xreg_coef,
                                                          xreg, n_cond);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // update variance estimate for the arima model - this was passed by reference
  sigma2 = exp(2 * solution.value);
  // do the same for model coefficients
  for (int i = 0; i < solution.x.size(); i++) {
    coef[i] = solution.x[i];
  }
  std::vector<double> solution_to_expand(coef.size());
  for (int i = 0; i < coef.size(); i++) {
    solution_to_expand[i] = coef[i];
  }
  // before we can initialize the state space model we need to transform the
  // parameters one last time
  arima_transform_parameters<seasonal, false>(solution_to_expand, kind);
  // finally, make state space model
  structural_model<double> arima_ss = make_arima( solution_to_expand,
                                                  delta, kind, kappa, ss_init);
  model.set(arima_ss);
}

#endif
