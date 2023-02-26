#ifndef ARIMA_CSS_SOLVER
#define ARIMA_CSS_SOLVER

#include "utils/utils.h"

#include "arima/structures/structural_model.h"
#include "arima/structures/fitting_method.h"
#include "arima/solvers/arima_gradient.h"

#include "arima/solvers/state_space.h"
#include "utils/xreg.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"

using FunctionXd = cppoptlib::function::Function<double>;

template <const bool has_xreg, const bool seasonal>
class ARIMA_CSS_PROBLEM : public FunctionXd {

  using EigVec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<double, Eigen::VectorXd, Eigen::MatrixXd>;

private:
  ArimaLossGradient<seasonal, has_xreg, double> Grad;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // initialize with a given arima structure
    ARIMA_CSS_PROBLEM(const std::vector<double> &y,
                      const arima_kind &kind,
                      const bool intercept,
                      const bool drift,
                      std::vector<std::vector<double>> &xreg) {
      this->Grad = ArimaLossGradient<seasonal, has_xreg>(
        y, kind, intercept,
        vec_to_mat(xreg, y.size(), intercept, drift), 1);
    }
    double operator()(const EigVec &x) {
      return this->Grad.loss(x);
    }
    StateXd Eval(const EigVec &x,
                 const int order = 1) {
      StateXd state(x.size(), 1);
      state.x = x;
      state.gradient = this->Grad.Gradient(x);
      state.value = this->Grad.loss(x);
      return state;
    }
    void finalize( structural_model<double> &model,
                   const EigVec & final_pars,
                   const double kappa,
                   const SSinit ss_init) {
      this->Grad.finalize(model, final_pars, kappa, ss_init);
    }
};

template <const bool has_xreg, const bool seasonal>
double arima_solver_css(std::vector<double> &y,
                      const arima_kind &kind,
                      structural_model<double> &model,
                      std::vector<std::vector<double>> xreg,
                      const bool intercept,
                      const bool drift,
                      std::vector<double> &coef,
                      const double kappa,
                      const SSinit ss_init) {

  Eigen::VectorXd x(coef.size());
  for(size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal>> solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal> css_arima_problem(
      y, kind, intercept, drift, xreg);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // create state space representation for model
  css_arima_problem.finalize( model, solution.x, kappa, ss_init);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < coef.size(); i++) coef[i] = solution.x[i];
  // return sigma2
  return exp(2 * solution.value);
}

#endif
