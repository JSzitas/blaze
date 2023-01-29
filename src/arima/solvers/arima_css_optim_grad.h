#ifndef ARIMA_CSS_SOLVER
#define ARIMA_CSS_SOLVER

#include "utils/utils.h"

#include "arima/structures/structural_model.h"
#include "arima/structures/fitting_method.h"
#include "arima/solvers/arima_gradient.h"

#include "arima/solvers/arima_css_likelihood.h"
#include "arima/solvers/state_space.h"

#include "arima/utils/xreg.h"

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
  arima_kind kind;
  ArimaLossGradient<seasonal, has_xreg, double> Grad;
  StateXd state;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    // initialize with a given arima structure
    ARIMA_CSS_PROBLEM(std::vector<double> &y,
                      const arima_kind &kind,
                      lm_coef<double> &xreg_pars,
                      std::vector<std::vector<double>> &xreg,
                      size_t n_cond) {
      this->kind = kind;
      // initialize an xreg matrix
      size_t n = y.size();
      std::vector<double> _xreg = flatten_vec(xreg);
      EigMat new_mat = Eigen::Map<EigMat>( _xreg.data(), n, xreg.size());
      if (xreg_pars.has_intercept()) {
        new_mat.conservativeResize(Eigen::NoChange, new_mat.cols() + 1);
        new_mat.col(new_mat.cols() - 1) = EigVec::Constant(n, 1, 1);
      }
      this->Grad = ArimaLossGradient<seasonal, has_xreg>(
        y, kind, xreg_pars.has_intercept(), new_mat, n_cond, 1);
      this->state = StateXd(kind.p() + kind.P() + kind.q() + kind.Q() +
        new_mat.cols(), 1);
    }
    double operator()(const EigVec &x) {
      return this->Grad.loss(x);
    }
    StateXd Eval(const Eigen::VectorXd &x,
                 const int order = 1) {

      this->state.x = x;
      this->state.gradient = this->Grad.Gradient(x);
      this->state.value = this->Grad.loss(x);
      return this->state;
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
  for (auto &val : x) {
    val = 0;
  }
  // initialize to all zeroes except for xreg
  for (size_t i = arma_size; i < vec_size; i++) {
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
  // std::cout << " Function evaluations taken: " << css_arima_problem.f_evals << std::endl;
  // std::cout << " Solver iterations: " << solver_state.num_iterations << std::endl;
  // solver_state.Status == IterationLimit ||      GradientNormViolation,  // Minimum norm in gradient vector has been reached.
  //     HessianConditionViolation  // Maximum condition number of hessian_t has been reached.
  // };

  // enum class Status {
  //   NotStarted = -1,
  //     Continue = 0,     // Optimization should continue.
  //     IterationLimit,   // Maximum of allowed iterations has been reached.
  //     XDeltaViolation,  // Minimum change in parameter vector has been reached.
  //     FDeltaViolation,  // Minimum chnage in cost function has been reached.
  //     GradientNormViolation,  // Minimum norm in gradient vector has been reached.
  //     HessianConditionViolation  // Maximum condition number of hessian_t has been
  //   // reached.
  // };
  // if we are to return standard errors as well
  // if constexpr(return_hessian) {
  //   // first get the computed numerical hessian
  //   auto state = css_arima_problem.Eval(solution.x);
  //   //std::cout << *(state.hessian) << std::endl;
  //   // available under state.hessian
  //   // next multiply by -1 and the number of available observations
  //   auto est_hessian = state.hessian * -1 * (double)n_available;
  //   // next, invert this
  //   // the result is the variance covariance matrix of coefficient estimates, except
  //   // for the fact that xreg coefficients (including intercept) have significantly
  //   // understated effects
  //   // since for typical use (e.g. coefficient standard errors) you only want the
  //   // diagonal elements anyways, you can simply multiply the diagonal entries
  //   // on those elements by the sample variance used in originally scaling the data
  // }
  css_arima_problem.finalize( model, solution.x, delta, kappa, ss_init);
  // pass fitted coefficients back to the caller
  for (size_t i = 0; i < vec_size; i++) {
    coef[i] = solution.x[i];
  }
}

#endif
