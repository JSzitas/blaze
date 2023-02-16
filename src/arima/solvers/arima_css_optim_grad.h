#ifndef ARIMA_CSS_SOLVER
#define ARIMA_CSS_SOLVER

#include "utils/utils.h"

#include "arima/structures/structural_model.h"
#include "arima/structures/fitting_method.h"
#include "arima/solvers/arima_gradient.h"

// #include "arima/solvers/arima_css_likelihood.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/xreg.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"
// profiler
// #include "third_party/coz.h"

using FunctionXd = cppoptlib::function::Function<double>;

template <const bool has_xreg, const bool seasonal>
class ARIMA_CSS_PROBLEM : public FunctionXd {

  using EigVec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<double, Eigen::VectorXd, Eigen::MatrixXd>;
private:
  const arima_kind kind;
  ArimaLossGradient<seasonal, has_xreg, double> Grad;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // initialize with a given arima structure
    ARIMA_CSS_PROBLEM(const std::vector<double> &y,
                      const arima_kind &kind,
                      const bool intercept,
                      std::vector<std::vector<double>> &xreg,
                      size_t n_cond) : kind(kind) {
      // initialize an xreg matrix
      size_t n = y.size();
      this->Grad = ArimaLossGradient<seasonal, has_xreg>(
        y, kind, intercept,
        vec_to_mat(xreg, y.size(), intercept), n_cond, 1);
    }
    double operator()(const EigVec &x) {
      return this->Grad.loss(x);
      // COZ_PROGRESS_NAMED("Function evaluation");
    }
    StateXd Eval(const Eigen::VectorXd &x,
                 const int order = 1) {
      StateXd state(x.size(), 1);
      state.x = x;
      state.gradient = this->Grad.Gradient(x);
      // COZ_PROGRESS_NAMED("Gradient evaluation");
      state.value = this->Grad.loss(x);
      // COZ_PROGRESS_NAMED("State evaluation");
      return state;
    }
    void finalize( structural_model<double> &model,
                   const Eigen::VectorXd & final_pars,
                   std::vector<double> & delta,
                   double kappa,
                   SSinit ss_init) {
      // this function creates state space representation of the ARIMA model
      auto x_temp = this->Grad.get_expanded_coef(final_pars);
      structural_model<double> arima_ss = make_arima( x_temp,
                                                      delta, this->kind,
                                                      kappa, ss_init);
      model.set(arima_ss);
      auto y_temp = this->Grad.get_y_temp(final_pars);
      // get arima steady state values
      arima_steady_state(y_temp, model);
    }
};

template <const bool has_xreg, const bool seasonal>
void arima_solver_css(std::vector<double> &y, structural_model<double> &model,
                      lm_coef<double> xreg_coef,
                      std::vector<std::vector<double>> xreg,
                      const arima_kind &kind, std::vector<double> &coef,
                      std::vector<double> &delta, const int n_cond,
                      const int n_available, const double kappa,
                      const SSinit ss_init, double &sigma2) {

  size_t vec_size = kind.p() + kind.q() + kind.P() + kind.Q() + xreg_coef.size();
  size_t arma_size = kind.p() + kind.q() + kind.P() + kind.Q();
  Eigen::VectorXd x(vec_size);
  for(size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // using xreg
  for (size_t i = arma_size; i < vec_size; i++)
    x[i] = xreg_coef.coef[i - arma_size];
  // initialize solver
  cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal>> solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal> css_arima_problem(
      y, kind, xreg_coef.has_intercept(), xreg, n_cond);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(css_arima_problem, x);
  // update variance estimate for the arima model - this was passed by reference
  sigma2 = exp(2 * solution.value);
  css_arima_problem.finalize( model, solution.x, delta, kappa, ss_init);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < vec_size; i++) {
    coef[i] = solution.x[i];
  }
}

#endif
