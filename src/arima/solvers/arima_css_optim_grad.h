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


template <const bool has_xreg, const bool seasonal, typename scalar_t=float>
class ARIMA_CSS_PROBLEM : public cppoptlib::function::Function<scalar_t, ARIMA_CSS_PROBLEM<has_xreg, seasonal,scalar_t>> {

  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<scalar_t, EigVec, EigMat>;

private:
  ArimaLossGradient<seasonal, has_xreg, scalar_t> Grad;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // initialize with a given arima structure
    ARIMA_CSS_PROBLEM(const std::vector<scalar_t> &y,
                      const arima_kind &kind,
                      const bool intercept,
                      const bool drift,
                      std::vector<std::vector<scalar_t>> &xreg) {
      this->Grad = ArimaLossGradient<seasonal, has_xreg, scalar_t>(
        y, kind, intercept,
        vec_to_mat(xreg, y.size(), intercept, drift), 1);
    }
    scalar_t operator()(const EigVec &x) {
      return this->Grad.loss(x);
    }
    // add impl of grad, hessian, eval
    void Gradient(const EigVec &x, EigVec *grad) {
      (*grad) = this->Grad.Gradient(x);
    }
    StateXd Eval(const EigVec &x,
                 const int order = 1) {
      StateXd state(x.size(), 1);
      state.x = x;
      state.gradient = this->Grad.Gradient(x);
      state.value = this->Grad.loss(x);
      return state;
    }
    void finalize( structural_model<scalar_t> &model,
                   const EigVec & final_pars,
                   const scalar_t kappa,
                   const SSinit ss_init) {
      this->Grad.finalize(model, final_pars, kappa, ss_init);
    }
};

template <const bool has_xreg, const bool seasonal, typename scalar_t=float>
scalar_t arima_solver_css(std::vector<scalar_t> &y,
                      const arima_kind &kind,
                      structural_model<scalar_t> &model,
                      std::vector<std::vector<scalar_t>> xreg,
                      const bool intercept,
                      const bool drift,
                      std::vector<scalar_t> &coef,
                      const scalar_t kappa,
                      const SSinit ss_init) {

  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> x(coef.size());
  for(size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<ARIMA_CSS_PROBLEM<has_xreg, seasonal, scalar_t>> solver;
  // and arima problem
  ARIMA_CSS_PROBLEM<has_xreg, seasonal, scalar_t> css_arima_problem(
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
