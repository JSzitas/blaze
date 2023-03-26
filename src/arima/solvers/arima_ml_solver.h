#ifndef ARIMA_ML_SOLVER
#define ARIMA_ML_SOLVER

#include "arima/solvers/arima_ml_kalman.h"

#include "utils/xreg.h"

// include optimizer library
#include "third_party/eigen.h"
#include "third_party/optim.h"


template <const SSinit ss_type, const bool has_xreg, const bool seasonal, const bool transform, typename scalar_t=float>
class ARIMA_ML_PROBLEM : public cppoptlib::function::Function<scalar_t,
                                                     ARIMA_ML_PROBLEM<ss_type, has_xreg, seasonal, transform, scalar_t>> {
  using EigVec = Eigen::Matrix<scalar_t, Eigen::Dynamic, 1>;
  using EigMat = Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic>;
  using StateXd = cppoptlib::function::State<scalar_t, EigVec, EigMat>;

  KalmanARIMA<ss_type, seasonal, has_xreg, transform, scalar_t> MLSolver;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_ML_PROBLEM<ss_type, has_xreg, seasonal, transform, scalar_t>(
      const std::vector<scalar_t> &y,
      const arima_kind &kind,
      const bool intercept,
      const bool drift,
      std::vector<std::vector<scalar_t>> &xreg,
      const scalar_t kappa) {
      this->MLSolver = KalmanARIMA<ss_type, seasonal, has_xreg,
        transform, scalar_t
      >(y, kind, vec_to_mat(xreg, y.size(), intercept, drift), kappa);
  }
  scalar_t operator()(const EigVec &x) {
    // print_vector(x);
   return this->MLSolver(x);
  }
  // add impl of grad, hessian, eval
  void Gradient(const EigVec &x, EigVec *grad) {
    cppoptlib::utils::ComputeFiniteGradient(*this, x, grad);
  }
  void Hessian(const EigVec &x, EigMat *hessian) {
    cppoptlib::utils::ComputeFiniteHessian(*this, x, hessian);
  }
  StateXd Eval(const EigVec &x,
               const int order = 1) {
    StateXd state(x.size(), 1);
    state.value = this->operator()(x);
    state.x = x;
    if (order >= 1) {
      this->Gradient(x, &state.gradient);
    }
    if ((order >= 2) && (state.hessian)) {
      this->Hessian(x, &*(state.hessian));
    }
    return state;
  }
  structural_model<scalar_t> get_structural_model() const {
    return this->MLSolver.get_structural_model();
  }
  scalar_t get_sigma() const {
    return this->MLSolver.get_sigma();
  }
};

template <const SSinit ss_type, const bool has_xreg, const bool seasonal,
          const bool transform, typename scalar_t=float>
scalar_t solver_ml(std::vector<scalar_t> &y,
                 structural_model<scalar_t> &model,
                 const bool intercept,
                 const bool drift,
                 std::vector<std::vector<scalar_t>> &xreg,
                 const arima_kind &kind,
                 std::vector<scalar_t> &coef,
                 const scalar_t kappa) {
  Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> x(coef.size());
  for (size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<
    ARIMA_ML_PROBLEM<ss_type,has_xreg, seasonal, transform, scalar_t>> solver;
  // and arima problem
  ARIMA_ML_PROBLEM<ss_type, has_xreg, seasonal, transform, scalar_t> ml_arima_problem(
      y, kind, intercept, drift, xreg, kappa);
  // and finally, minimize
  auto [solution, solver_state] = solver.Minimize(ml_arima_problem, x);
  // replace model V
  // ml_arima_problem.finalize();
  // write back structural model
  model = ml_arima_problem.get_structural_model();
  for (size_t i = 0; i < coef.size(); i++) coef[i] = solution.x[i];
  // return variance estimate for the arima model
  return ml_arima_problem.get_sigma();
}


template <const bool has_xreg, const bool seasonal, const bool transform,
          typename scalar_t=float>
scalar_t arima_solver_ml(std::vector<scalar_t> &y,
                       structural_model<scalar_t> &model,
                       const bool intercept,
                       const bool drift,
                       std::vector<std::vector<scalar_t>> &xreg,
                       const arima_kind &kind,
                       std::vector<scalar_t> &coef,
                       const scalar_t kappa,
                       const SSinit ss_init) {
  // dispatch on SSinit
  if( ss_init == SSinit::Gardner ) {
    return solver_ml<SSinit::Gardner, has_xreg, seasonal, transform, scalar_t>(
        y, model, intercept, drift, xreg, kind, coef, kappa);
  }
  else if( ss_init == SSinit::Rossignol) {
    return solver_ml<SSinit::Rossignol, has_xreg, seasonal, transform, scalar_t>(
        y, model, intercept, drift, xreg, kind, coef, kappa);
  }
  return 0.0;
}

#endif
