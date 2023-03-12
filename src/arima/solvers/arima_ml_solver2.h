#ifndef ARIMA_ML_SOLVER
#define ARIMA_ML_SOLVER

#include "arima/solvers/arima_ml_kalman.h"

#include "third_party/eigen.h"
#include "third_party/optim.h"

#include "utils/xreg.h"

using FunctionXd = cppoptlib::function::Function<double>;

template <const SSinit ss_type, const bool has_xreg, const bool seasonal, const bool transform>
class ARIMA_ML_PROBLEM : public FunctionXd {
  KalmanARIMA<ss_type, seasonal, has_xreg, transform, true, double> MLSolver;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // initialize with a given arima structure
  ARIMA_ML_PROBLEM<ss_type, has_xreg, seasonal, transform>(
      const std::vector<double> &y,
      const arima_kind &kind,
      const bool intercept,
      const bool drift,
      std::vector<std::vector<double>> &xreg,
      const double kappa) {
      this->MLSolver = KalmanARIMA<ss_type, seasonal, has_xreg, transform
        >(y, kind, vec_to_mat(xreg, y.size(), intercept, drift), kappa);
  }
  double operator()(const Eigen::VectorXd &x) {
   return this->MLSolver(x);
  }
  structural_model<double> get_structural_model() const {
    return this->MLSolver.get_structural_model();
  }
  double get_sigma() const {
    return this->MLSolver.get_sigma();
  }
};

template <const SSinit ss_type, const bool has_xreg, const bool seasonal, const bool transform>
double solver_ml(std::vector<double> &y,
                 structural_model<double> &model,
                 const bool intercept,
                 const bool drift,
                 std::vector<std::vector<double>> &xreg,
                 const arima_kind &kind,
                 std::vector<double> &coef,
                 const double kappa) {
  Eigen::VectorXd x(coef.size());
  for (size_t i = 0; i < coef.size(); i++) x(i) = coef[i];
  // initialize solver
  cppoptlib::solver::Bfgs<
    ARIMA_ML_PROBLEM<ss_type,has_xreg, seasonal, transform>> solver;
  // and arima problem
  ARIMA_ML_PROBLEM<ss_type, has_xreg, seasonal, transform> ml_arima_problem(
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


template <const bool has_xreg, const bool seasonal, const bool transform>
double arima_solver_ml(std::vector<double> &y,
                       structural_model<double> &model,
                       const bool intercept,
                       const bool drift,
                       std::vector<std::vector<double>> &xreg,
                       const arima_kind &kind,
                       std::vector<double> &coef,
                       const double kappa,
                       const SSinit ss_init) {
  // dispatch on SSinit
  if( ss_init == SSinit::Gardner ) {
    return solver_ml<SSinit::Gardner, has_xreg, seasonal, transform>(
        y, model, intercept, drift, xreg, kind, coef, kappa);
  }
  else if( ss_init == SSinit::Rossignol) {
    return solver_ml<SSinit::Rossignol, has_xreg, seasonal, transform>(
        y, model, intercept, drift, xreg, kind, coef, kappa);
  }
  return 0.0;
}

#endif
